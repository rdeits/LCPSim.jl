using RigidBodyDynamics
using MeshCat
using MeshCatMechanisms
using LCPSim
using JuMP, Gurobi
using Polyhedra, CDDLib
using Base.Test

function box_demo()
    urdf = joinpath(@__DIR__, "box.urdf")
    urdf_mech = parse_urdf(Float64, urdf)
    mechanism, base = planar_revolute_base()
    attach!(mechanism, base, urdf_mech)
    world = root_body(mechanism)

    vis = Visualizer()
    if !haskey(ENV, "CI")
        open(vis)
        wait(vis)
    end
    mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf), vis[:robot])
    floor = planar_obstacle(default_frame(world), [0, 0, 1.], [0, 0, 0.])
    walls = [
        planar_obstacle(default_frame(world), [1., 0, 0], [-1., 0, 0]),
        planar_obstacle(default_frame(world), [-1., 0, 0], [1., 0, 0]),
    ]

    boundary = SimpleHRepresentation(vcat(eye(3), -eye(3)), vcat([1.1, 1.1, 2.0], -[-1.1, -1.1, -0.1]))
    for (i, obstacle) in enumerate([floor, walls...])
        A = vcat([h.outward_normal.v' for h in obstacle.interior]...)
        b = vcat([h.outward_normal.v' * h.point.v for h in obstacle.interior]...)
        interior = SimpleHRepresentation(A, b)
        p = intersect(boundary, interior)
        setobject!(vis[:environment][string(i)], CDDPolyhedron{3, Float64}(p))
    end

    core = findbody(mechanism, "core")

    env = Environment([
        (core, pt, obs) for pt in [
            Point3D(default_frame(core), 0.1, 0, 0.2),
            Point3D(default_frame(core), -0.1, 0, 0.2),
            Point3D(default_frame(core), 0.1, 0, -0.2),
            Point3D(default_frame(core), -0.1, 0, -0.2)
        ] for obs in [floor, walls...]
    ])

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

        results = LCPSim.simulate(x0, controller, env, Δt, N, GurobiSolver(Gurobi.Env(), OutputFlag=0));
        for r in results
            set_configuration!(mvis, configuration(r.state))
            sleep(Δt * 0.01)
        end
        @test length(results) == N
    end
end

box_demo()
