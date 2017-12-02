using DrakeVisualizer
DrakeVisualizer.any_open_windows() || DrakeVisualizer.new_window()
using RigidBodyDynamics
using RigidBodyTreeInspector
using LCPSim
using JuMP, Gurobi
using Polyhedra, CDDLib
using Base.Test

function DrakeVisualizer.addgeometry!(vis::Visualizer, obs::Obstacle, boundary::HRepresentation)
    p = intersect(boundary, obs.interior)
    addgeometry!(vis, CDDPolyhedron{3, Float64}(p))
end

function box_demo()
    urdf_mech = parse_urdf(Float64, joinpath(@__DIR__, "box.urdf"))
    mechanism, base = planar_revolute_base()
    attach!(mechanism, base, urdf_mech)
    world = root_body(mechanism)

    basevis = Visualizer()[:box]
    delete!(basevis)
    vis = basevis[:robot]
    setgeometry!(vis, mechanism, parse_urdf(joinpath(@__DIR__, "box.urdf"), mechanism))

    floor = planar_obstacle(default_frame(world), [0, 0, 1.], [0, 0, 0.])
    walls = [
        planar_obstacle(default_frame(world), [1., 0, 0], [-1., 0, 0]),
        planar_obstacle(default_frame(world), [-1., 0, 0], [1., 0, 0]),
    ]

    boundary = SimpleHRepresentation(vcat(eye(3), -eye(3)), vcat([1.1, 1.1, 2.0], -[-1.1, -1.1, -0.1]))
    addgeometry!.(basevis[:environment], [floor, walls...], boundary)

    core = findbody(mechanism, "core")

    env = Environment(
        Dict(core => ContactEnvironment(
                    [
                    Point3D(default_frame(core), 0.1, 0, 0.2),
                    Point3D(default_frame(core), -0.1, 0, 0.2),
                    Point3D(default_frame(core), 0.1, 0, -0.2),
                    Point3D(default_frame(core), -0.1, 0, -0.2),
                     ],
                    [floor, walls...])));

    controller = x -> zeros(num_velocities(x))
    Δt = 0.02
    N = 100

    for i in 1:50
        @show i
        srand(i)
        x0 = MechanismState{Float64}(mechanism)
        set_velocity!(x0, zeros(num_velocities(x0)))
        set_configuration!(x0, findjoint(mechanism, "base_x"), rand(1) .- 0.5)
        set_configuration!(x0, findjoint(mechanism, "base_z"), [1.0])
        set_configuration!(x0, findjoint(mechanism, "base_rotation"), randn(1))
        set_velocity!(x0, findjoint(mechanism, "base_x"), 2 * randn(1))
        set_velocity!(x0, findjoint(mechanism, "base_z"), 2 * randn(1))
        
        results = LCPSim.simulate(x0, controller, env, Δt, N, GurobiSolver(OutputFlag=0));
        for r in results
            set_configuration!(x0, configuration(r.state))
            settransform!(vis, x0)
            sleep(Δt * 0.01)
        end
        @test length(results) == N
    end
end

box_demo()
