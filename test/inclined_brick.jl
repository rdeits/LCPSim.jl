using LCPSim
using Base.Test
using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using StaticArrays: SVector
# using Cbc: CbcSolver
using Gurobi: GurobiSolver
using Rotations: RotY

urdf = joinpath(@__DIR__, "..", "examples", "box.urdf")

function inclined_brick(θ)
    mechanism = parse_urdf(Float64, urdf)
    core = findbody(mechanism, "core")
    fixed_joint = joint_to_parent(core, mechanism)
    floating_base = Joint(fixed_joint.name, frame_before(fixed_joint), frame_after(fixed_joint), 
                          Planar([1., 0, 0], [0., 0, 1.]),
                          position_bounds=[Bounds(-5., 5), Bounds(0., 3), Bounds(-2π, 2π)],
                          velocity_bounds=[Bounds(-10., 10), Bounds(-10., 10), Bounds(-2π, 2π)],
                          effort_bounds=[Bounds(0., 0) for i in 1:3])
    replace_joint!(mechanism, fixed_joint, floating_base)

    world = root_body(mechanism)
    R = RotY(θ)
    floor = planar_obstacle(default_frame(world), R * SVector(0., 0, 1), [0, 0, 0.], 1.0)
    env = Environment(
        Dict(core => ContactEnvironment(
                    [
                    Point3D(default_frame(core), SVector(0.1, 0, 0.2)),
                    Point3D(default_frame(core), SVector(-0.1, 0, 0.2)),
                    Point3D(default_frame(core), SVector(0.1, 0, -0.2)),
                    Point3D(default_frame(core), SVector(-0.1, 0, -0.2)),
                     ],
                    [floor]
                    )))

    x0 = MechanismState{Float64}(mechanism)
    set_velocity!(x0, zeros(num_velocities(x0)))

    center = R * (SVector(-0.8, 0, 0.1))
    set_configuration!(x0, floating_base, [center[1], center[3], -θ + π/2])
    mechanism, env, x0
end


@testset "bricks on inclined planes" begin
    # A brick on an incline which is slightly too shallow, so the brick does not
    # slide at all
    mechanism, env, x1 = inclined_brick(π/4 - 0.1)
    Δt = 0.05
    controller = x -> zeros(num_velocities(x))
    q0 = copy(configuration(x1))
    results_stick = LCPSim.simulate(x1, controller, env, Δt, 50, GurobiSolver(OutputFlag=0))
    @test length(results_stick) == 50 
    for i in 8:length(results_stick)
        @test norm(configuration(results_stick[i].state) .- configuration(results_stick[i - 1].state)) <= 1e-5
        @test norm(velocity(results_stick[i].state)) <= 1e-5
    end

    # A slightly steeper incline causes the brick to begin to slide
    mechanism, env, x2 = inclined_brick(π/4 + 0.1)
    Δt = 0.05
    controller = x -> zeros(num_velocities(x))
    q0 = copy(configuration(x2))
    results_slide = LCPSim.simulate(x2, controller, env, Δt, 50, GurobiSolver(OutputFlag=0))
    @test length(results_slide) == 50 
    for i in 8:length(results_slide)
        @test norm(configuration(results_slide[i].state) .- configuration(results_slide[i - 1].state)) >= 1e-2
        @test norm(velocity(results_slide[i].state)) >= 0.5
    end

    if Pkg.installed("RigidBodyTreeInspector") !== nothing
        @eval using RigidBodyTreeInspector
        @eval using DrakeVisualizer; 
        @eval using CoordinateTransformations
        @eval using GeometryTypes
        DrakeVisualizer.any_open_windows() || DrakeVisualizer.new_window()

        v1 = Visualizer()[:box][:stick]
        setgeometry!(v1, mechanism, parse_urdf(urdf, mechanism))
        v2 = Visualizer()[:box][:slide]
        setgeometry!(v2, mechanism, parse_urdf(urdf, mechanism))

        for i in 1:length(results_stick)
            set_configuration!(x1, configuration(results_stick[i].state))
            settransform!(v1, x1)
            set_configuration!(x2, configuration(results_slide[i].state))
            settransform!(v2, x2)
            sleep(Δt)
        end
    end

end
