module ContactLQR

using RigidBodyDynamics
using ForwardDiff

"""
`care(A, B, Q, R)`

Compute 'X', the solution to the continuous-time algebraic Riccati equation,
defined as A'X + XA - (XB)R^-1(B'X) + Q = 0, where R is non-singular.

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf

Implementation from https://github.com/JuliaControl/ControlSystems.jl/blob/master/src/matrix_comps.jl
by Jim Crist and other contributors.
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

Implementation from https://github.com/JuliaControl/ControlSystems.jl/blob/master/src/synthesis.jl
by Jim Crist and other contributors.
"""
function lqr(A, B, Q, R)
    S = care(A, B, Q, R)
    K = R\B'*S
    return K, S
end

"""
    `dare(A, B, Q, R)`

Compute `X`, the solution to the discrete-time algebraic Riccati equation,
defined as A'XA - X - (A'XB)(B'XB + R)^-1(B'XA) + Q = 0, where A and R
are non-singular.

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf

Implementation from https://github.com/JuliaControl/ControlSystems.jl/blob/master/src/matrix_comps.jl
by Jim Crist and other contributors.
"""
function dare(A, B, Q, R)
    G = try
        B*inv(R)*B'
    catch
        error("R must be non-singular.")
    end

    Ait = try
        inv(A)'
    catch
        error("A must be non-singular.")
    end

    Z = [A + G*Ait*Q   -G*Ait;
         -Ait*Q        Ait]

    S = schurfact(Z)
    S = ordschur(S, abs.(S.values).<=1)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]
    return U21/U11
end

"""
    `dlqr(A, B, Q, R)`

Calculate the optimal gain matrix `K` for the state-feedback law `u[k] = K*x[k]` that
minimizes the cost function:

J = sum(x'Qx + u'Ru, 0, inf).

For the discrete time model `x[k+1] = Ax[k] + Bu[k]`.

Implementation from https://github.com/JuliaControl/ControlSystems.jl/blob/master/src/synthesis.jl
by Jim Crist and other contributors.
```
"""
function dlqr(A, B, Q, R)
    S = dare(A, B, Q, R)
    K = (B'*S*B + R)\(B'S*A)
    return K, S
end


function contact_jacobian(state, contacts)
    if isempty(contacts)
        return zeros(0, num_positions(state))
    end
    q = configuration(state)
    contact_jacobians = Matrix{eltype(q)}[]
    for contact in contacts
        J = ForwardDiff.jacobian(q) do q
            x = MechanismState{eltype(q)}(state.mechanism)
            set_configuration!(x, q)
            T = transform_to_root(x, contact.frame)
            (T * contact).v
        end
        push!(contact_jacobians, J)
    end
    Jc = vcat(contact_jacobians...)
    Jc = Jc[[i for i in 1:size(Jc, 1) if !all(Jc[i, :] .== 0)], :]
    Jc
end

function dynamics_with_contact_constraint(state::MechanismState, input::AbstractVector, Jc::AbstractMatrix)
    v = velocity(state)
    M = full(mass_matrix(state))
    C_plus_g = dynamics_bias(state)
    St = eye(num_velocities(state))
    for joint in joints(state.mechanism)
        bounds = effort_bounds(joint)
        for (i, b) in enumerate(bounds)
            if b.upper == b.lower == 0
                j = parentindexes(velocity(state, joint))[i]
                St[j, j] = 0
            end
        end
    end
    Jct = Jc'
    Jcbar = inv(M) * Jct * inv(Jc * inv(M) * Jct)
    Φ = (I - Jcbar * Jc) * inv(M)
    ϕ = -Φ * (C_plus_g)

    v̇ = Φ * St * input + ϕ
    vcat(v, v̇)
end

"""
From "Full Dynamics LQR Control of a Humanoid Robot: An Experimental Study on
Balancing and Squatting" by Sean Mason et al.
"""
function contact_linearize(state0, input0, Jc)
    if norm(velocity(state0)) > 0
        error("Only static postures supported")
    end
    mechanism = state0.mechanism

    function dynamics(x, u)
        q = x[1:num_positions(state0)]
        v = x[num_positions(state0) + 1 : end]
        state = MechanismState(mechanism, q, v)
        dynamics_with_contact_constraint(state, u, Jc)
    end

    A = ForwardDiff.jacobian(state_vector(state0)) do x
        dynamics(x, input0)
    end
    B = ForwardDiff.jacobian(input0) do u
        dynamics(state_vector(state0), u)
    end
    A, B, dynamics(state_vector(state0), input0)
end

"""
    A_d, B_d, c_d = zero_order_hold(A, B, c, Δt)

Discretize the continuous-time linear system:

    ẋ = Ax + Bu + c

into the discrete-time linear system:

    x[i+1] = Ax[i] + Bu[i] + c

using a zero-order (piecewise constant) hold on u(t) with discretization
interval Δt.
"""
function zero_order_hold(A, B, c, Δt)
    nx = size(A, 1)
    nu = size(B, 2)
    matc = zeros(nx + nu + 1, nx + nu + 1)
    matc[1:nx, :] = hcat(A, B, c)
    matd = expm(matc * Δt)
    A_d = matd[1:nx, 1:nx]
    B_d = matd[1:nx, nx + (1:nu)]
    c_d = matd[1:nx, nx + nu + 1]
    A_d, B_d, c_d
end

"""
From "Balancing and Walking Using Full Dynamics LQR Control With Contact
Constraints" by Sean Mason et al.
"""
function contact_lqr(state::MechanismState, input::AbstractVector, Q::AbstractMatrix, R::AbstractMatrix, contacts::AbstractVector{<:Point3D})
    Jc = contact_jacobian(state, contacts)
    A, B, c = contact_linearize(state, input, Jc)
    A[abs.(A) .< 1e-6] .= 0
    N = nullspace([Jc zeros(Jc); zeros(Jc) Jc])
    Am = N' * A * N
    Bm = N' * B
    Rm = R
    Qm = N' * Q * N
    Km, Sm = lqr(Am, Bm, Qm, Rm)
    K = Km * N'
    S = N * Sm * N'
    return K, S
end

function contact_dlqr(state::MechanismState, input::AbstractVector, Q::AbstractMatrix, R::AbstractMatrix, Δt, contacts::AbstractVector{<:Point3D}=[])
    Jc = contact_jacobian(state, contacts)
    A, B, c = contact_linearize(state, input, Jc)
    A[abs.(A) .< 1e-6] .= 0
    N = nullspace([Jc zeros(Jc); zeros(Jc) Jc])
    Am = N' * A * N
    Bm = N' * B
    Rm = R
    Qm = N' * Q * N
    Am_d, Bm_d, cm_d = zero_order_hold(Am, Bm, zeros(size(Am, 1)), Δt)
    Km_d, Sm_d = dlqr(Am_d, Bm_d, Qm, Rm)
    K_d = Km_d * N'
    S_d = N * Sm_d * N'
    return K_d, S_d
end

end
