using LCPSim
using Base.Test
using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using StaticArrays: SVector
using Gurobi: GurobiSolver
using Rotations: RotY
using MeshCat
using MeshCatMechanisms
using CoordinateTransformations: Translation

function load_model(urdf, θ)
    mechanism = parse_urdf(Float64, urdf)
    env = LCPSim.parse_contacts(mechanism, urdf, 1.0, :xz)
    planar_joint = findjoint(mechanism, "plane_to_box")
    planar_joint.position_bounds .= Bounds(-10, 10)
    planar_joint.velocity_bounds .= Bounds(-1000, 1000)
    state = MechanismState(mechanism)
    set_configuration!(state, planar_joint, [-1, -1.5, θ])
    mechanism, env, state
end

@testset "bricks on inclined planes" begin
    vis = Visualizer()
    if !haskey(ENV, "CI")
        open(vis)
        wait(vis)
    end

    # A brick on an incline which is slightly too shallow, so the brick does not
    # slide at all
    urdf = joinpath(@__DIR__, "urdf", "shallow_incline.urdf")
    mechanism, env, x1 = load_model(urdf, -π/4 - 0.1)
    mv1 = MechanismVisualizer(mechanism, URDFVisuals(urdf), vis[:stick][:model])
    Δt = 0.02
    N = 100
    controller = x -> zeros(num_velocities(x))

    q0 = copy(configuration(x1))
    results_stick = LCPSim.simulate(x1, controller, env, Δt, N, GurobiSolver(Gurobi.Env(), OutputFlag=0))
    @testset "sticking" begin
        @test length(results_stick) == N
        for i in 50:length(results_stick)
            @test norm(configuration(results_stick[i].state) .- configuration(results_stick[i - 1].state)) <= 1e-5
            @test norm(velocity(results_stick[i].state)) <= 1e-5
        end
    end

    # A slightly steeper incline causes the brick to begin to slide
    urdf = joinpath(@__DIR__, "urdf", "steep_incline.urdf")
    mechanism, env, x2 = load_model(urdf, -π/4 + 0.1)
    mv2 = MechanismVisualizer(mechanism, URDFVisuals(urdf), vis[:slide][:model])
    settransform!(vis[:slide], Translation(0, 1, 0))
    q0 = copy(configuration(x2))
    results_slide = LCPSim.simulate(x2, controller, env, Δt, N, GurobiSolver(Gurobi.Env(), OutputFlag=0))
    @testset "sliding" begin
        @test length(results_slide) == N
        for i in 50:length(results_slide)
            @test norm(configuration(results_slide[i].state) .- configuration(results_slide[i - 1].state)) >= 1e-2
            @test norm(velocity(results_slide[i].state)) >= 0.5
        end
    end

    for i in 1:length(results_stick)
        set_configuration!(mv1, configuration(results_stick[i].state))
        set_configuration!(mv2, configuration(results_slide[i].state))
        sleep(Δt)
    end

end
