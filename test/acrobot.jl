using Base.Test
using RigidBodyDynamics
using RigidBodyDynamics: Bounds
using StaticArrays: SVector
using Cbc: CbcSolver
using LCPSim
using MeshCat
using MeshCatMechanisms


#=
"""
`care(A, B, Q, R)`

Compute 'X', the solution to the continuous-time algebraic Riccati equation,
defined as A'X + XA - (XB)R^-1(B'X) + Q = 0, where R is non-singular.

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf
"""
function care(A, B, Q, R)
    G = try
        B*inv(R)*B'
    catch
        error("R must be non-singular.")
    end

    Z = [A  -G;
        -Q  -A']

    S = schurfact(Z)
    S = ordschur(S, real(S.values).<0)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]
    return U21/U11
end

"""
`lqr(A, B, Q, R)`

Calculate the optimal gain matrix `K` for the state-feedback law `u = K*x` that
minimizes the cost function:

J = integral(x'Qx + u'Ru, 0, inf).

For the continuous time model `dx = Ax + Bu`.

`lqr(sys, Q, R)`

Solve the LQR problem for state-space system `sys`. Works for both discrete
and continuous time systems.

See also `LQG`

Usage example:
```julia
A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = eye(2)
R = eye(1)
L = lqr(sys,Q,R)

u(t,x) = -L*x # Form control law,
t=0:0.1:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0)
plot(t,x, lab=["Position", "Velocity"]', xlabel="Time [s]")
```
"""
function lqr(A, B, Q, R)
    S = care(A, B, Q, R)
    K = R\B'*S
    return K
end

function dynamics(mechanism, x, u)
    state = MechanismState{eltype(x)}(mechanism)
    result = DynamicsResult{eltype(x)}(mechanism)
    ẋ = zeros(eltype(x), 4)
    dynamics!(ẋ, result, state, x, u)
    ẋ
end

A = ForwardDiff.jacobian(x -> dynamics(mechanism, x, zeros(eltype(x), 2)), [π, 0, 0, 0])
B = ForwardDiff.jacobian(u -> dynamics(mechanism, eltype(u)[π, 0, 0, 0], [0, u[2]]), [0., 0])
Q = diagm([10, 10, 1, 1])
R = diagm([1, 1])
K = lqr(A, B, Q, R)
=#

@testset "acrobot" begin
    K = [0.0 0.0 0.0 0.0; -266.087 -107.349 -114.95 -54.3175]
    urdf = joinpath(@__DIR__, "..", "examples", "urdf", "Acrobot.urdf")
    mechanism = parse_urdf(Float64, urdf)
    world = root_body(mechanism)
    env = Environment{Float64}([])


    x0 = MechanismState{Float64}(mechanism)
    set_velocity!(x0, zeros(num_velocities(x0)))
    set_configuration!(x0, findjoint(mechanism, "shoulder"), [π])
    set_configuration!(x0, findjoint(mechanism, "elbow"), [-0.1])
    q0 = copy(configuration(x0))
    v0 = copy(velocity(x0))

    controller = x -> begin
        x̂ = vcat(configuration(x), velocity(x)) - vcat(q0, v0)
        -K * x̂
    end

    Δt = 0.05
    N = 100
    results = LCPSim.simulate(x0, controller, env, Δt, N, CbcSolver());
    @test length(results) == N
    @test norm(velocity(results[end].state)) < 0.03
    @test norm(configuration(results[end].state) - q0) < 1.2


    vis = Visualizer()
    if !haskey(ENV, "CI")
        open(vis)
        wait(vis)
    end

    mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf), vis)
    for r in results
        sleep(Δt)
        set_configuration!(mvis, configuration(r.state))
    end
end
