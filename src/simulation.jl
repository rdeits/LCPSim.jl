

function add_free_region_constraints!(model::Model, xnext::LinearizedState, env::Environment)
    for (body, contact_env) in env.contacts
        for contact_point in contact_env.points
            # function _point_in_world(x)
            #     q = x[1:num_positions(xnext)]
            #     v = x[(num_positions(xnext)+1):end]
            #     x_diff = MechanismState(xnext.mechanism, q, v)
            #     point_in_world = transform_to_root(x_diff, contact_point.frame) * contact_point
            #     point_in_world.v
            # end

            # position = Point3D(root_frame(xnext.mechanism), _point_in_world(state_vector(x_dynamics)) + ForwardDiff.jacobian(_point_in_world, state_vector(x_dynamics)) * (state_vector(xnext) - state_vector(x_dynamics)))

            position_in_world = Linear.evaluate(x -> transform_to_root(x, contact_point.frame) * contact_point, xnext)

            ConditionalJuMP.disjunction!(model,
                [@?(position_in_world ∈ P) for P in contact_env.free_regions]) # (7)
        end
    end
end


function update(x::MechanismState{X, M}, 
                u, 
                joint_limits::Associative{<:Joint, <:HRepresentation}, 
                env::Environment, 
                Δt::Real, 
                model::Model, 
                x_dynamics::MechanismState{<:Number}=x) where {X, M}
    mechanism = x.mechanism
    world = root_body(mechanism)
    qnext = @variable(model, [1:num_positions(x)], lowerbound=-10, basename="qnext", upperbound=10)
    vnext = @variable(model, [1:num_velocities(x)], lowerbound=-10, basename="vnext", upperbound=10)
    xnext = LinearizedState(x_dynamics, vcat(qnext, vnext))
    # xnext = MechanismState(mechanism, qnext, vnext)

    contact_results = map(env.contacts) do item
        body, contact_env = item
        body => [resolve_contact(xnext, body, contact_point, obstacle, model)
            for contact_point in contact_env.points for obstacle in contact_env.obstacles]
    end

    externalwrenches = map(contact_results) do item
        body, results = item
        body => sum([Wrench(transform_to_root(x_dynamics, result.point.frame) * result.point,
                    contact_force(result)) for result in results])
    end

    add_free_region_constraints!(model, xnext, env)

    # function _config_derivative(v)
    #     q = oftype(v, configuration(x_dynamics))
    #     x_diff = MechanismState(mechanism, q, v)
    #     configuration_derivative(x_diff)
    # end
    # jac_dq_wrt_v = ForwardDiff.jacobian(_config_derivative, velocity(x_dynamics))

    jac_dq_wrt_v = Linear.jacobian(configuration_derivative, xnext)[:, length(qnext) + 1:end]



    joint_limit_results = convert(Dict{Joint, Vector{JointLimitResult{Variable, M}}},
        Dict([joint => resolve_joint_limits(xnext, joint, limits, model) for (joint, limits) in joint_limits]))
    # joint_limit_results = Dict{Joint, Vector{JointLimitResult{Variable, Vector{GenericAffExpr{M, Variable}}}}}([joint => resolve_joint_limits(xnext, joint, limits, model) for (joint, limits) in joint_limits])
    joint_limit_forces = zeros(GenericAffExpr{M, Variable}, num_velocities(x))
    for (joint, results) in joint_limit_results
        for result in results
            joint_limit_forces .+= (jac_dq_wrt_v')[:, parentindexes(configuration(x, joint))...] * generalized_force(result)
        end
    end

    H = mass_matrix(x_dynamics)
    bias = dynamics_bias(x_dynamics, externalwrenches)
    config_derivative = jac_dq_wrt_v * vnext

    @constraint(model, H * (vnext - velocity(x)) .== Δt * (u .+ joint_limit_forces .- bias)) # (5)
    @constraint(model, qnext .- configuration(x) .== Δt .* config_derivative) # (6)

    LCPUpdate(MechanismState(mechanism, qnext, vnext), u, contact_results, joint_limit_results)
end

function simulate(x0::MechanismState, 
                  controller, 
                  joint_limits::Associative{<:Joint, <:HRepresentation}, 
                  env::Environment, 
                  Δt::Real, 
                  N::Integer,
                  solver::JuMP.MathProgBase.SolverInterface.AbstractMathProgSolver)
    x = x0
    map(1:N) do i
        m = Model(solver=solver)
        u = controller(x)
        up = update(x, u, joint_limits, env, Δt, m)
        solve(m)
        update_value = getvalue(up)
        x = update_value.state
        update_value
    end
end

function optimize(x0::MechanismState, 
                  input_limits::Associative{<:Joint, <:AbstractVector},
                  joint_limits::Associative{<:Joint, <:HRepresentation}, 
                  env::Environment, 
                  Δt,
                  N::Integer,
                  m::Model=Model())
    x = x0
    u_min = zeros(num_velocities(x0))
    u_max = zeros(num_velocities(x0))
    for (joint, limits) in input_limits
        u_min[parentindexes(velocity(x0, joint))...] .= limits[1]
        u_max[parentindexes(velocity(x0, joint))...] .= limits[2]
    end
    results = map(1:N) do i
        u = @variable(m, [1:num_velocities(x0)], basename="u_$i")
        setlowerbound.(u, u_min)
        setupperbound.(u, u_max)
        for i in 1:num_velocities(x0)
            if u_min[i] == u_max[i]
                JuMP.fix(u[i], u_min[i])
            end
        end
        up = update(x, u, joint_limits, env, Δt, m, x0)
        x = up.state
        up
    end
    m, results
end
