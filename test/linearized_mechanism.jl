using MechanismComplementarity
using RigidBodyDynamics
using Base.Test

const lineval = MechanismComplementarity.Linear.evaluate

srand(1)
@testset "linearized mechanism" begin
    @testset "truly linear mechanism" begin
        mechanism = rand_chain_mechanism(Float64, [Prismatic{Float64} for i = 1:5]...)
        for i in 1:10
            x0 = MechanismState(mechanism, randn(num_positions(mechanism)), randn(num_velocities(mechanism)))
            x = MechanismState(mechanism, randn(num_positions(mechanism)), randn(num_velocities(mechanism)))
            x_linear = MechanismComplementarity.Linear.LinearizedState(x0, state_vector(x))
            @test center_of_mass(x) ≈ lineval(center_of_mass, x_linear)

            for body in bodies(mechanism)
                frame = default_frame(body)
                @test transform_to_root(x, frame) ≈ lineval(x -> transform_to_root(x, frame), x_linear)

                p = Point3D(frame, randn(3))
                @test transform_to_root(x, frame) * p ≈ lineval(x -> transform_to_root(x, frame) * p, x_linear)
            end
        end
    end
end


    # world = RigidBody{Float64}("world")
    # brick = Mechanism(world; gravity=SVector(0, 0, -9.81))

    # frame = CartesianFrame3D("dummy")
    # inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    # dummy = RigidBody(inertia)
    # base_x = Joint("base_x", Prismatic([1., 0, 0]))
    # attach!(brick, world, base_x, eye(Transform3D, frame_before(base_x), default_frame(world)), dummy)

    # frame = CartesianFrame3D("core")
    # inertia = SpatialInertia(frame, 0.1 * eye(3), zeros(3), 1.0)
    # core = RigidBody(inertia)
    # base_z = Joint("base_z", Prismatic([0., 0, 1]))
    # attach!(brick, dummy, base_z, eye(Transform3D, frame_before(base_z), default_frame(dummy)), core)