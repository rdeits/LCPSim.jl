using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using MeshCat
using MeshCatMechanisms
using LCPSim
using JuMP, Gurobi
using Base.Test

@testset "box in a box" begin
    urdf = joinpath(@__DIR__, "urdf", "box_in_a_box.urdf")
    mechanism = parse_urdf(Float64, urdf)
    planar_joint = findjoint(mechanism, "floor_to_box")
    planar_joint.position_bounds .= [
        Bounds(-3, 3),
        Bounds(-3, 3),
        Bounds(-2π, 2π)
    ]
    planar_joint.velocity_bounds .= Bounds(-1000, 1000)

    vis = Visualizer()
    if !haskey(ENV, "CI")
        open(vis)
        wait(vis)
    end
    mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf), vis[:robot])
    env = LCPSim.parse_contacts(mechanism, urdf, 1.0, :xz)
    controller = x -> zeros(num_velocities(x))
    Δt = 0.02
    N = 100
    solver = GurobiSolver(Gurobi.Env(), OutputFlag=0)

    for i in 1:50
        @show i
        srand(i)
        x0 = MechanismState{Float64}(mechanism)
        set_velocity!(x0, zeros(num_velocities(x0)))
        set_configuration!(x0, planar_joint, [rand() - 0.5, -1.0, randn()])
        set_velocity!(x0, planar_joint, 2 * randn(3))
        results = LCPSim.simulate(x0, controller, env, Δt, N, solver);
        for r in results
            set_configuration!(mvis, configuration(r.state))
            sleep(0.1 * Δt)
        end
        @test length(results) == N
    end
end
