module MechanismComplementarity

using Polyhedra
using StaticArrays
using JuMP
using JuMP: GenericAffExpr
using ConditionalJuMP
using Base.Test
using RigidBodyDynamics
using Rotations
using ForwardDiff

module Linear

    using RigidBodyDynamics
    using ForwardDiff

    struct LinearizedState{S <: MechanismState, S2 <: MechanismState, SDual <: MechanismState{<:ForwardDiff.Dual}}
        linear_state::S
        current_state::S2
        dual_state::SDual
    end

    function LinearizedState(linear::MechanismState, current::MechanismState)
        nq = num_positions(linear)
        nv = num_velocities(linear)
        na = num_additional_states(linear)
        N = nq + nv + na
        xdiff = Vector{ForwardDiff.Dual{N, Float64}}(N)
        ForwardDiff.seed!(xdiff, state_vector(linear), ForwardDiff.construct_seeds(ForwardDiff.Partials{N, Float64}))
        diffstate = MechanismState(linear.mechanism, xdiff[1:nq], xdiff[nq + (1:nv)], xdiff[nq + nv + (1:na)])
        LinearizedState(linear, current, diffstate)
    end

    unwrap(p::Point3D) = (v -> Point3D(p.frame, v), p.v)
    unwrap(p::FreeVector3D) = (v -> FreeVector3D(p.frame, v), p.v)
    unwrap(p) = (identity, p)

    function evaluate(f::Function, s::LinearizedState)
        wrapper, ydual = unwrap(f(s.dual_state))
        nx = length(state_vector(s.current_state))
        v = ForwardDiff.value.(ydual)
        
        if isa(v, AbstractVector)
            J = similar(v, (length(v), nx))
            ForwardDiff.extract_jacobian!(J, ydual, nx)
            wrapper(v .+ J * (state_vector(s.current_state) - state_vector(s.linear_state)))
        else
            wrapper(v + ForwardDiff.partials(ydual)' * (state_vector(s.current_state) - state_vector(s.linear_state)))
        end 
    end

    function jacobian(f::Function, s::LinearizedState)
        wrapper, ydual = unwrap(f(s.dual_state))
        nx = length(state_vector(s.current_state))
        v = ForwardDiff.value.(ydual)
        J = similar(v, (length(v), nx))
        ForwardDiff.extract_jacobian!(J, ydual, nx)
        J
    end
        


end

using .Linear: LinearizedState

include("environments.jl")
include("conditional_jump_extensions.jl")
include("variable_containers.jl")
include("contact.jl")
include("joint_limits.jl")
include("simulation.jl")

end # module
