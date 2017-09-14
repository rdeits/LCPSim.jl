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
    for (body, contact_env) in env.contacts
        contact_results[body] = []
        for obstacle in contact_env.obstacles, contact_point in contact_env.points
            push!(contact_results[body], resolve_contact(xnext, body, contact_point, obstacle, model))
        end
    end

    externalwrenches = map(contact_results) do item
        body, results = item
        body => sum([Wrench(transform_to_root(x_dynamics, result.point.frame) * result.point,
                    contact_force(result)) for result in results])
    end

    jac_dq_wrt_v = Linear.jacobian(configuration_derivative, xnext)[:, length(qnext) + 1:end]

    joint_limit_results::Dict{Joint, Vector{JointLimitResult{Variable, M}}} = 
        Dict([joint => resolve_joint_limits(xnext, joint, model) for joint in joints(mechanism)])

    joint_limit_forces = zeros(GenericAffExpr{M, Variable}, num_velocities(x))
    for (joint, results) in joint_limit_results
        for result in results
            joint_limit_forces .+= (jac_dq_wrt_v')[:, parentindexes(configuration(xnext.dual_state, joint))...] * generalized_force(result)
        end
    end

    H = mass_matrix(x_dynamics)
    bias = dynamics_bias(x_dynamics, externalwrenches)
    config_derivative = jac_dq_wrt_v * vnext

    @constraint(model, H * (vnext - velocity(x)) .== Δt .* (u .+ joint_limit_forces .- bias)) # (5)
    @constraint(model, qnext .- configuration(x) .== Δt .* config_derivative) # (6)

    LCPUpdate(StateRecord(mechanism, vcat(qnext, vnext)), u, contact_results, joint_limit_results)
end

"""
Update the linearization state of ``xnext`` using the current state ``x``
"""
function semi_implicit_update!(xnext::LinearizedState, x::StateRecord, Δt)
    set_linearization_velocity!(xnext, velocity(x))
    qnext = configuration(xnext.linearization_state) .+ Δt .* configuration_derivative(xnext.linearization_state)
    set_linearization_configuration!(xnext, qnext)
end

function simulate(x0::MechanismState{T, M}, 
                  controller, 
                  env::Environment, 
                  Δt::Real, 
                  N::Integer,
                  solver::JuMP.MathProgBase.SolverInterface.AbstractMathProgSolver) where {T, M}
    x = StateRecord(x0)
    xnext = LinearizedState{Variable}(x0)
    input_limits = all_effort_bounds(x0.mechanism)
    results = LCPUpdate{Float64, M, Vector{Float64}}[]
    for i in 1:N
        m = Model(solver=solver)
        u = clamp.(controller(x), input_limits)
        semi_implicit_update!(xnext, x, Δt)
        # xnext.linearization_state.q[1:4] = normalize(xnext.linearization_state.q[1:4])
        # x.configuration[1:4] = normalize(x.configuration[1:4])
        up = update(x, xnext, u, env, Δt, m)
        status = solve(m)
        if status != :Optimal
            break
        end
        update_value = getvalue(up)
        x = update_value.state
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
        set_linearization_configuration!(xnext, configuration(seed[i].state))
        set_linearization_velocity!(xnext, velocity(seed[i].state))
        u = @variable(m, [1:num_velocities(x0)], basename="u_$i")
        setbounds.(u, input_limits)
        fix_if_tightly_bounded.(u)
        up = update(x, xnext, u, env, Δt, m)
        x = up.state
        up
    end
    setvalue.(results, seed)
    m, results
end
    
