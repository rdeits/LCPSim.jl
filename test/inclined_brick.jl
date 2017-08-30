using LCPSim
using Base.Test
using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using StaticArrays: SVector
using Cbc: CbcSolver
using Rotations: RotY

const urdf = joinpath(@__DIR__, "..", "examples", "box.urdf")

function inclined_brick(θ)
    mechanism = parse_urdf(Float64, urdf)
    core = findbody(mechanism, "core")
    floating_base = joint_to_parent(core, mechanism)
    floating_base.jointType = Planar([1., 0, 0], [0., 0, 1.])
    floating_base.position_bounds = [Bounds(-5., 5), Bounds(0., 3), Bounds(-2π, 2π)]
    floating_base.velocity_bounds = [Bounds(-10., 10), Bounds(-10., 10), Bounds(-2π, 2π)]
    floating_base.effort_bounds = [Bounds(0., 0) for i in 1:3]

    world = root_body(mechanism)
    R = RotY(θ)
    floor = planar_obstacle(default_frame(world), R * SVector(0., 0, 1), [0, 0, 0.], 1.0)
    free_space = space_between([floor])
    env = Environment(
        Dict(core => ContactEnvironment(
                    [
                    Point3D(default_frame(core), SVector(0.1, 0, 0.2)),
                    Point3D(default_frame(core), SVector(-0.1, 0, 0.2)),
                    Point3D(default_frame(core), SVector(0.1, 0, -0.2)),
                    Point3D(default_frame(core), SVector(-0.1, 0, -0.2)),
                     ],
                    [floor],
                    [free_space])))

    x0 = MechanismState{Float64}(mechanism)
    set_velocity!(x0, zeros(num_velocities(x0)))

    center = R * (SVector(-0.8, 0, 0.1))
    set_configuration!(x0, floating_base, [center[1], center[3], -θ + π/2])
    mechanism, env, x0
end


# A brick on an incline which is slightly too shallow, so the brick does not
# slide at all
@testset "sticking brick" begin
    mechanism, env, x0 = inclined_brick(π/4 - 0.01)
    Δt = 0.05
    controller = x -> zeros(num_velocities(x))
    q0 = copy(configuration(x0))
    results = LCPSim.simulate(x0, controller, env, Δt, 15, CbcSolver())
    for r in results
        @test configuration(r.state) ≈ q0
        @test isapprox(velocity(r.state), [0., 0, 0], atol=1e-11)
    end
end

# A slightly steeper incline causes the brick to begin to slide
@testset "sliding brick" begin
    mechanism, env, x0 = inclined_brick(π/4 + 0.01)
    Δt = 0.05
    controller = x -> zeros(num_velocities(x))
    q0 = copy(configuration(x0))
    results = LCPSim.simulate(x0, controller, env, Δt, 15, CbcSolver())
    @test !(configuration(results[end].state) ≈ q0)
    @test velocity(results[end].state)[2] < -0.1
end


