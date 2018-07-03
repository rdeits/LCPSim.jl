using RigidBodyDynamics: Bounds, upper, lower

all_configuration_bounds(m::Mechanism) =
    collect(Base.Iterators.flatten(map(position_bounds, joints(m))))

all_velocity_bounds(m::Mechanism) =
    collect(Base.Iterators.flatten(map(velocity_bounds, joints(m))))

all_effort_bounds(m::Mechanism) =
    collect(Base.Iterators.flatten(map(effort_bounds, joints(m))))

function setbounds(x::Variable, b::Bounds)
    setlowerbound(x, lower(b))
    setupperbound(x, upper(b))
    x
end

function create_state_variables(model::Model, mech::Mechanism)
    qnext::Vector{Variable} = @variable(model, [1:num_positions(mech)], basename="qnext")
    vnext::Vector{Variable} = @variable(model, [1:num_velocities(mech)], basename="vnext")

    setbounds.(qnext, all_configuration_bounds(mech))
    setbounds.(vnext, all_velocity_bounds(mech))

    @assert all(isfinite, lowerbound.(qnext)) && all(isfinite, upperbound.(qnext)) "All joint position limits must be finite"
    @assert all(isfinite, lowerbound.(vnext)) && all(isfinite, upperbound.(vnext)) "All joint velocity limits must be finite"
    qnext, vnext
end

function update(x::StateRecord{X, M},
                xnext::LinearizedState{Variable},
                u,
                env::Environment,
                Δt::Real,
                model::Model) where {X, M}
    x_dynamics = linearization_state(xnext)
    mechanism = x.mechanism
    qnext, vnext = create_state_variables(model, mechanism)
    set_current_configuration!(xnext, qnext)
    set_current_velocity!(xnext, vnext)

    contact_results = Dict{RigidBody{M}, Vector{ContactResult{Variable, M}}}()
    for (body, contact_point, obstacle) in env.contacts
        results = get!(() -> [], contact_results, body)
        push!(results, resolve_contact(xnext, body, contact_point, obstacle, model))
    end

    world = root_body(mechanism)
    externalwrenches = Dict{RigidBody{M}, Wrench{GenericAffExpr{M, Variable}}}()
    for (body, results) in contact_results
        for result in results
            c = transform(x_dynamics, contact_force(result), default_frame(world))
            p = transform(x_dynamics, result.point, default_frame(world))
            w = Wrench(p, c)
            if haskey(externalwrenches, body)
                externalwrenches[body] += w
            else
                externalwrenches[body] = w
            end
        end
    end

    joint_limit_results = [resolve_joint_limit(model, xnext, joint) for joint in joints(mechanism)]
    joint_limit_forces = zeros(GenericAffExpr{M, Variable}, num_velocities(x))

    for (i, joint) in enumerate(joints(mechanism))
        vrange = velocity_range(x_dynamics, joint)
        jt = joint_type(joint)
        jac_v_wrt_dq = RigidBodyDynamics.configuration_derivative_to_velocity_jacobian(jt, configuration(x_dynamics, joint))
        joint_limit_forces[velocity_range(x_dynamics, joint)] .+= jac_v_wrt_dq * generalized_force(joint_limit_results[i])
    end

    bias_coriolis_gravity = linearized(dynamics_bias, xnext)

    τ_external_wrenches = zeros(GenericAffExpr{M, Variable}, num_velocities(x_dynamics))
    for (body, wrench) in externalwrenches
        J = geometric_jacobian(x_dynamics, path(mechanism, body, world))
        τ_external_wrenches .+= torque(J, wrench)
    end

    bias = bias_coriolis_gravity + τ_external_wrenches
    Q_v = RigidBodyDynamics.velocity_to_configuration_derivative_jacobian(x_dynamics)
    config_derivative = Q_v * vnext

    H = mass_matrix(x_dynamics)
    HΔv = H * (vnext - velocity(x))
    @constraint(model, HΔv .== Δt .* (u .+ joint_limit_forces .- bias)) # (5)
    @constraint(model, qnext .- configuration(x) .== Δt .* config_derivative) # (6)

    LCPUpdate(StateRecord(mechanism, vcat(qnext, vnext)), u, contact_results, joint_limit_results)
end

"""
Update the linearization state of ``xnext`` using the current state ``x``
"""
function semi_implicit_update!(xnext::LinearizedState, x::Union{StateRecord, MechanismState}, Δt)
    # we need to use xnext to compute qdot instead of x since
    # x might be a StateRecord and not a full MechanismState
    set_linearization_configuration!(xnext, configuration(x))
    set_linearization_velocity!(xnext, velocity(x))
    q̇ = configuration_derivative(linearization_state(xnext))

    set_linearization_configuration!(xnext,
        configuration(x) .+ Δt .* q̇)
end

function explicit_update!(xnext::LinearizedState, x::StateRecord)
    set_linearization_velocity!(xnext, velocity(x))
    set_linearization_configuration!(xnext, configuration(x))
end


function simulate(x0::MechanismState{T, M},
                  controller,
                  env::Environment,
                  Δt::Real,
                  N::Integer,
                  solver::JuMP.MathProgBase.SolverInterface.AbstractMathProgSolver;
                  termination::Function = state -> false,
                  relinearize=true) where {T, M}
    x = StateRecord(x0)
    xnext = LinearizedState{Variable}(x0)
    input_limits = all_effort_bounds(x0.mechanism)
    results = LCPUpdate{Float64, M, Float64}[]
    for i in 1:N
        m = Model(solver=solver)
        if relinearize
            semi_implicit_update!(xnext, x, Δt)
            # explicit_update!(xnext, x)
        end
        u = clamp.(controller(x), input_limits)
        up = update(x, xnext, u, env, Δt, m)
        if i > 1
            setvalue(up, results[end])
            ConditionalJuMP.warmstart!(m, false)
        end
        status = solve(m; suppress_warnings=true)
        if status != :Optimal
            break
        end
        update_value = getvalue(up)
        x = update_value.state
        if termination(x)
            break
        end
        push!(results, update_value)
    end
    results
end

function fix_if_tightly_bounded(x::Variable)
    if getlowerbound(x) == getupperbound(x)
        JuMP.fix(x, getlowerbound(x))
    end
end

function optimize(x0::MechanismState,
                  env::Environment,
                  Δt,
                  N::Integer,
                  m::Model=Model())
    x = StateRecord(x0)
    xnext = LinearizedState{Variable}(x0)
    input_limits = all_effort_bounds(x0.mechanism)
    results = map(1:N) do i
        u = @variable(m, [1:num_velocities(x0)], basename="u_$i")
        setbounds.(u, input_limits)
        fix_if_tightly_bounded.(u)
        up = update(x, xnext, u, env, Δt, m)
        x = up.state
        up
    end
    m, results
end

function optimize(x0::MechanismState,
                  env::Environment,
                  Δt,
                  seed::AbstractVector{<:LCPUpdate},
                  m::Model=Model())
    N = length(seed)
    x = StateRecord(x0)
    xnext = LinearizedState{Variable}(x0)
    input_limits = all_effort_bounds(x0.mechanism)
    results = map(1:N) do i
        if i > 1
            semi_implicit_update!(xnext, seed[i - 1].state, Δt)
        else
            semi_implicit_update!(xnext, x, Δt)
        end
        u = @variable(m, [1:num_velocities(x0)], basename="u_$i")
        setbounds.(u, input_limits)
        fix_if_tightly_bounded.(u)
        up = update(x, xnext, u, env, Δt, m)
        x = up.state
        up
    end
    setvalue.(results, seed)
    ConditionalJuMP.warmstart!(m, false)
    m, results
end

