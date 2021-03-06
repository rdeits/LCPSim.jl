# Comparison of two (equivalent) ways of constructing a mechanism with an
# x, y, θ base: an explicit Planar joint vs. two Prismatic joints and one Rotational
# joint in series.


using LCPSim
using Base.Test
using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using StaticArrays: SVector
using Cbc: CbcSolver
using Gurobi: GurobiSolver
using Rotations: RotY
using MeshCat
using MeshCatMechanisms

urdf = joinpath(@__DIR__, "urdf", "box.urdf")

function box_with_planar_base()
    mechanism = parse_urdf(Float64, urdf)
    core = findbody(mechanism, "core")
    fixed_joint = joint_to_parent(core, mechanism)
    floating_base = Joint(fixed_joint.name, frame_before(fixed_joint), frame_after(fixed_joint),
                          Planar([1., 0, 0], [0., 0, 1.]),
                          position_bounds=[Bounds(-5., 5), Bounds(0., 3), Bounds(-2π, 2π)],
                          velocity_bounds=[Bounds(-100., 10), Bounds(-100., 10), Bounds(-100., 100)],
                          effort_bounds=[Bounds(0., 0) for i in 1:3])
    replace_joint!(mechanism, fixed_joint, floating_base)

    x0 = MechanismState{Float64}(mechanism)
    set_velocity!(x0, zeros(num_velocities(x0)))
    set_configuration!(x0, floating_base, [-1, 0.3, -0.3])
    configuration_derivative_to_velocity!(velocity(x0, floating_base), floating_base, configuration(x0, floating_base), [3., 0, 0])
    mechanism, x0
end

function box_with_dummy_links()
    urdf_mech = parse_urdf(Float64, urdf)
    mechanism, base = planar_revolute_base()
    attach!(mechanism, base, urdf_mech)

    x0 = MechanismState{Float64}(mechanism)
    set_velocity!(x0, zeros(num_velocities(x0)))
    set_configuration!(x0, findjoint(mechanism, "base_x"), [-1])
    set_configuration!(x0, findjoint(mechanism, "base_z"), [0.3])
    set_configuration!(x0, findjoint(mechanism, "base_rotation"), [0.3])
    set_velocity!(x0, findjoint(mechanism, "base_x"), [3])
    mechanism, x0
end

function environment_with_floor(mechanism)
    world = root_body(mechanism)
    floor = planar_obstacle(default_frame(world), [0, 0, 1.], [0, 0, 0.], 0.5)

    core = findbody(mechanism, "core")
    env = Environment([
        (core, pt, floor) for pt in [
                    Point3D(default_frame(core), SVector(0.1, 0, 0.2)),
                    Point3D(default_frame(core), SVector(-0.1, 0, 0.2)),
                    Point3D(default_frame(core), SVector(0.1, 0, -0.2)),
                    Point3D(default_frame(core), SVector(-0.1, 0, -0.2)),
                     ]
            ])
    env
end

function withenv(f::Function)
    env = Gurobi.Env()
    try
        f(env)
    finally
        # https://github.com/JuliaOpt/Gurobi.jl/issues/110
        # Gurobi.free_env(env)
    end
end

@testset "planar vs dummy link brick" begin
    mech1, x1 = box_with_planar_base()
    env1 = environment_with_floor(mech1)

    mech2, x2 = box_with_dummy_links()
    env2 = environment_with_floor(mech2)

    controller = LCPSim.passive_controller()
    Δt = 0.01
    N = 600

    results1 = withenv() do env
        LCPSim.simulate(x1, controller, env1, Δt, N, GurobiSolver(env, OutputFlag=0))
    end
    results2 = withenv() do env
        LCPSim.simulate(x2, controller, env2, Δt, N, GurobiSolver(env, OutputFlag=0))
    end

    vis = Visualizer()
    if !haskey(ENV, "CI")
        open(vis)
        wait(vis)
    end

    mv1 = MechanismVisualizer(mech1, URDFVisuals(urdf), vis[:planar])
    mv2 = MechanismVisualizer(mech2, URDFVisuals(urdf), vis[:dummy])

    for i in 1:length(results1)
        set_configuration!(mv1, configuration(results1[i].state))
        set_configuration!(mv2, configuration(results2[i].state))
        sleep(Δt)
    end

    @test length(results1) == N
    @test length(results2) == N

    for i in 1:N
        set_configuration!(x1, configuration(results1[i].state))
        set_velocity!(x1, velocity(results1[i].state))
        set_configuration!(x2, configuration(results2[i].state))
        set_velocity!(x2, velocity(results2[i].state))
        T1 = transform_to_root(x1, findbody(mech1, "core"))
        T2 = transform_to_root(x2, findbody(mech2, "core"))
        atol = N <= 30 ? 1e-2 : 5e-2
        @test rotation(T1) ≈ rotation(T2) atol=atol
        @test translation(T1) ≈ translation(T2) atol=atol
    end
end
