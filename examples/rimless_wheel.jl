using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using MeshCat
using MeshCatMechanisms
using LCPSim
using JuMP, Gurobi
using Base.Test

@testset "rimless wheel" begin
    urdf = joinpath(@__DIR__, "urdf", "rimless_wheel.urdf")
    mechanism = parse_urdf(Float64, urdf)
    floating_base = findjoint(mechanism, "floating_base")
    position_bounds(floating_base) .= Bounds(-100, 100)
    velocity_bounds(floating_base) .= Bounds(-100, 100)

    vis = Visualizer()
    if !haskey(ENV, "CI")
        open(vis)
        wait(vis)
    end
    mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf), vis[:robot])

    x0 = MechanismState(mechanism)
    set_configuration!(x0, [-1.5, 1.8, -0.05])
    env = LCPSim.parse_contacts(mechanism, urdf);
    controller = x -> zeros(num_velocities(x))
    Δt = 0.005
    N = 1700
    lcp_solver = GurobiSolver(Gurobi.Env(); OutputFlag=0)

    results = LCPSim.simulate(x0, controller, env, Δt, N, lcp_solver)

    for r in results
        set_configuration!(mvis, configuration(r.state))
        sleep(Δt)
    end
end
