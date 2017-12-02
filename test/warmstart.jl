using LCPSim
using Base.Test
using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using RigidBodyTreeInspector
using Cbc: CbcSolver
using Gurobi: GurobiSolver
using Rotations: RotY
using JuMP: setvalue, getvalue

urdf = joinpath(@__DIR__, "..", "examples", "box.urdf")

function box_with_planar_base()
    mechanism = parse_urdf(Float64, urdf)
    core = findbody(mechanism, "core")
    fixed_joint = joint_to_parent(core, mechanism)
    floating_base = Joint(fixed_joint.name, frame_before(fixed_joint), frame_after(fixed_joint), 
                          Planar([1., 0, 0], [0., 0, 1.]),
                          position_bounds=[Bounds(-5., 5), Bounds(0., 3), Bounds(-2π, 2π)],
                          velocity_bounds=[Bounds(-10., 10), Bounds(-10., 10), Bounds(-10., 10)],
                          effort_bounds=[Bounds(0., 0) for i in 1:3])
    replace_joint!(mechanism, fixed_joint, floating_base)
    x0 = MechanismState{Float64}(mechanism)
    set_configuration!(x0, floating_base, [-1, 0.3, -0.3])
    configuration_derivative_to_velocity!(velocity(x0, floating_base), floating_base, configuration(x0, floating_base), [1., 0, 0])
    mechanism, x0
end

function environment_with_floor(mechanism)
    world = root_body(mechanism)
    floor = planar_obstacle(default_frame(world), [0, 0, 1.], [0, 0, 0.], 0.5)

    core = findbody(mechanism, "core")
    env = Environment(
        Dict(core => ContactEnvironment(
                    [
                    Point3D(default_frame(core), 0.1, 0, 0.2),
                    Point3D(default_frame(core), -0.1, 0, 0.2),
                    Point3D(default_frame(core), 0.1, 0, -0.2),
                    Point3D(default_frame(core), -0.1, 0, -0.2),
                     ],
                    [floor],
                    )))
    env
end

@testset "simulation warmstarts" begin
    mech1, x1 = box_with_planar_base()
    q0 = copy(configuration(x1))
    v0 = copy(velocity(x1))
    env1 = environment_with_floor(mech1)

    controller = x -> zeros(num_velocities(x))
    Δt = 0.01
    N = 50
    eps = 2e-3

    basevis = Visualizer()[:box]
    delete!(basevis)
    vis = basevis[:robot]
    setgeometry!(vis, mech1, parse_urdf(urdf, mech1))

    @testset "linear dynamics" begin
        # First test the linear dynamics
        set_configuration!(x1, q0)
        set_velocity!(x1, v0)
        results1 = LCPSim.simulate(x1, controller, env1, Δt, N, GurobiSolver(Gurobi.Env(), OutputFlag=0); relinearize=false)

        for r in results1
            set_configuration!(x1, configuration(r.state))
            settransform!(vis, x1)
            sleep(Δt)
        end

        # Construct an optimization linearized around the initial state
        set_configuration!(x1, q0)
        set_velocity!(x1, v0)
        model2, results2 = LCPSim.optimize(x1, env1, Δt, N)
        setvalue.(results2, results1)
        ConditionalJuMP.warmstart!(model2, false)

        # Test that all the optimization constraints are satisfied by the
        # simulation result
        for c in model2.linconstr
            @test c.lb - eps <= getvalue(c.terms) <= c.ub + eps
        end
    end

    @testset "nonlinear dynamics" begin
        # Now simulate the nonlinear dynamics
        set_configuration!(x1, q0)
        set_velocity!(x1, v0)
        results3 = LCPSim.simulate(x1, controller, env1, Δt, N, GurobiSolver(Gurobi.Env(), OutputFlag=0); relinearize=true)

        for r in results3
            set_configuration!(x1, configuration(r.state))
            settransform!(vis, x1)
            sleep(Δt)
        end

        # Construct an optimization linearzed around the trajectory
        set_configuration!(x1, q0)
        set_velocity!(x1, v0)
        model4, results4 = LCPSim.optimize(x1, env1, Δt, results3)
        setvalue.(results4, results3)
        ConditionalJuMP.warmstart!(model4, false)

        # Test that all the optimization constraints are satisfied by the
        # simulation result
        for c in model4.linconstr
            @test c.lb - eps <= getvalue(c.terms) <= c.ub + eps
            if !(c.lb - eps <= getvalue(c.terms) <= c.ub + eps)
                @show c
            end
        end
    end
end