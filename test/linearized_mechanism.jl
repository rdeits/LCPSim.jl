using LCPSim
using RigidBodyDynamics
using Base.Test

srand(1)
@testset "linearized mechanism" begin
    @testset "truly linear mechanism" begin
        mechanism = rand_chain_mechanism(Float64, [Prismatic{Float64} for i = 1:5]...)
        x0 = MechanismState{Float64}(mechanism)
        x = MechanismState{Float64}(mechanism)
        for i in 1:10
            set_configuration!(x0, randn(num_positions(x0)))
            set_velocity!(x0, randn(num_velocities(x0)))
            set_configuration!(x, randn(num_positions(x0)))
            set_velocity!(x, randn(num_velocities(x0)))
            x_linear = LinearizedState(x0, x)
            @test center_of_mass(x) ≈ linearized(center_of_mass, x_linear)
            @test !(center_of_mass(x) ≈ center_of_mass(x0))

            for body in bodies(mechanism)
                frame = default_frame(body)
                @test transform_to_root(x, frame) ≈ linearized(x -> transform_to_root(x, frame), x_linear)

                p = Point3D(frame, randn(3))
                @test transform_to_root(x, frame) * p ≈ linearized(x -> transform_to_root(x, frame) * p, x_linear)
            end
        end
    end

    @testset "nonlinear mechanism" begin
        mechanism = rand_chain_mechanism(Float64, [QuaternionFloating{Float64}; [Revolute{Float64} for i = 1:5]]...)
        x0 = MechanismState{Float64}(mechanism)
        x = MechanismState{Float64}(mechanism)
        for i in 1:10
            set_configuration!(x0, randn(num_positions(x0)))
            set_velocity!(x0, randn(num_velocities(x0)))

            set_configuration!(x, configuration(x0) + 1e-5 * (rand(num_positions(x0)) .- 0.5))
            set_velocity!(x, velocity(x0) + 1e-5 * (rand(num_velocities(x0)) .- 0.5))

            x_linear =  LCPSim.Linear.LinearizedState(x0, x)

            @test center_of_mass(x) ≈ linearized(center_of_mass, x_linear)
            @test !(center_of_mass(x0) ≈ center_of_mass(x))  # make sure we've perturbed enough for linearization to help

            for (i, body) in enumerate(bodies(mechanism))
                frame = default_frame(body)
                @test transform_to_root(x, frame) ≈ linearized(x -> transform_to_root(x, frame), x_linear)

                if i > 1
                    @test !(transform_to_root(x0, frame) ≈ transform_to_root(x, frame))
                end

                p = Point3D(frame, randn(3))
                @test transform_to_root(x, frame) * p ≈ linearized(x -> transform_to_root(x, frame) * p, x_linear)
            end
        end
    end
end
