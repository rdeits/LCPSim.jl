using LCPSim
using Base.Test
using ForwardDiff
using RigidBodyDynamics
using RigidBodyDynamics: Bounds, position_bounds, velocity_bounds, effort_bounds
using StaticArrays: SVector

function hopper()
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world; gravity=SVector(0, 0, -9.81))

    frame = CartesianFrame3D("core")
    inertia = SpatialInertia(frame, 0.1 * eye(3), zeros(3), 1.0)
    core = RigidBody(inertia)
    base_z = Joint("base_z", Prismatic([0., 0, 1]))
    position_bounds(base_z) .= Bounds(-10., 10)
    velocity_bounds(base_z) .= Bounds(-10., 10)
    effort_bounds(base_z) .= Bounds(0., 0)
    attach!(mechanism, world, core, base_z)

    frame = CartesianFrame3D("foot")
    inertia = SpatialInertia(frame, 0.02 * eye(3), zeros(3), 0.2)
    foot = RigidBody(inertia)
    leg_z = Joint("leg_z", Prismatic([0., 0, 1]))
    position_bounds(leg_z) .= Bounds(-10., 10)
    velocity_bounds(leg_z) .= Bounds(-10., 10)
    effort_bounds(leg_z) .= Bounds(-20., 20)
    attach!(mechanism, core, foot, leg_z)
    return mechanism
end

function hopper_dynamics_in_contact(x, u)
    q = x[1:2]
    v = x[3:4]
    # assume the foot is stuck to the ground, then we have:
    v̇ = zeros(promote_type(eltype(x), eltype(u)), 2)
    v̇[1] = -u[1] - 9.81
    v̇[2] = -v̇[1]
    vcat(v, v̇)
end

@testset "contact LQR" begin
    x_fp = [1., -1, 0, 0]
    u_fp = [-9.81]
    @test hopper_dynamics_in_contact(x_fp, u_fp) ≈ zeros(4)

    A = ForwardDiff.jacobian(x -> hopper_dynamics_in_contact(x, u_fp), x_fp)
    B = ForwardDiff.jacobian(u -> hopper_dynamics_in_contact(x_fp, u), u_fp)
    Q = diagm([100, 10, 1, 1])
    R = 0.1 * eye(1)
    Jc = [1. 1]
    Jcdot = [0. 0]
    N = nullspace([Jc zeros(1, 2); Jcdot Jc])
    Am = N' * A * N
    Bm = N' * B
    Rm = R
    Qm = N' * Q * N
    _, _, _, Km, Sm = LCPSim.ContactLQR.lqr(Am, Bm, Qm, Rm)
    K = Km * N'
    @test K ≈ [-16.5831 16.5831 -4.64576 4.64576] atol=1e-4

    mechanism = hopper()
    state = MechanismState(mechanism, x_fp[1:2], x_fp[3:4])
    foot = findbody(mechanism, "foot")
    contacts = [Point3D(default_frame(foot), SVector(0., 0, 0))]
    _, _, _, K2, S2 = LCPSim.ContactLQR.contact_lqr(state, [0, -9.81], Q, 0.1 * eye(2), contacts)

    @test K2[2, :] ≈ K[:]
end

