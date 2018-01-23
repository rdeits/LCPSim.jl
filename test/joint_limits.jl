using LCPSim
using Base.Test
using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using StaticArrays: SVector
using Gurobi

@testset "joint limits" begin
    @testset "1D mechanism" begin
        world = RigidBody{Float64}("world")
        mech = Mechanism(world; gravity=SVector(0, 0, -9.81))
        frame = CartesianFrame3D("link")
        inertia = SpatialInertia(frame, eye(3), SVector(0., 0, 0), 1.0)
        link = RigidBody(inertia)
        slider = Joint("slider", Prismatic(SVector(0., 0, 1)),
                       position_bounds=[Bounds(0.0, 1.0)],
                       velocity_bounds=[Bounds(-10., 10.)],
                       effort_bounds=[Bounds(0., 0.)])
        attach!(mech, world, link, slider)
        env = Environment{Float64}(Dict())
        x0 = MechanismState{Float64}(mech)
        set_configuration!(x0, [0.5])

        Δt = 0.05
        N = 100
        results = LCPSim.simulate(x0, x -> zeros(length(velocity(x))),
                                  env, Δt, N, GurobiSolver(Gurobi.Env(), OutputFlag=0))
        @test length(results) == N
        for r in results
            @test -1e-3 <= configuration(r.state)[1] <= 1 + 1e-3
        end
        @test configuration(results[end].state)[1] <= 1e-3
        @test norm(velocity(results[end].state)) <= 1e-3
    end

    @testset "planar mechanism" begin
        world = RigidBody{Float64}("world")
        mech = Mechanism(world; gravity=SVector(0, 0, -9.81))
        frame = CartesianFrame3D("link")
        inertia = SpatialInertia(frame, eye(3), SVector(0., 0, 0), 1.0)
        link = RigidBody(inertia)

        planar = Joint("planar", Planar([1., 0, 0], [0., 0, 1.]),
                       position_bounds=[Bounds(-5., 5), Bounds(0., 3), Bounds(-2π, 2π)],
                       velocity_bounds=[Bounds(-10., 10), Bounds(-10., 10), Bounds(-2π, 2π)],
                       effort_bounds=[Bounds(0., 0) for i in 1:3])
        attach!(mech, world, link, planar)
        env = Environment{Float64}(Dict())
        x0 = MechanismState{Float64}(mech)
        set_configuration!(x0, [0.5, 1.0, 0.3])

        Δt = 0.05
        N = 100
        results = LCPSim.simulate(x0, x -> zeros(length(velocity(x))),
                                  env, Δt, N, GurobiSolver(Gurobi.Env(), OutputFlag=0))
        @test length(results) == N
        for r in results
            @test -5.0 - 1e-3 <= configuration(r.state)[1] <= 5.0 + 1e-3
            @test 0.0 - 1e-3 <= configuration(r.state)[2] <= 3.0 + 1e-3
            @test -2π - 1e-3 <= configuration(r.state)[3] <= 2π + 1e-3
        end
        @test configuration(results[end].state) ≈ [0.5, 0.0, 0.3]
        @test norm(velocity(results[end].state)) ≈ 0 atol=1e-8
    end
end
