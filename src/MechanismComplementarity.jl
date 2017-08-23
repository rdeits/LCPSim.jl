__precompile__()

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
    import RigidBodyDynamics: configuration, velocity
    using ForwardDiff

    struct LinearizedState{T, S <: MechanismState{<:Number}, SDual <: MechanismState{<:ForwardDiff.Dual}}
        linearization_state::S
        current_state_vector::Vector{T}
        dual_state::SDual
    end

    function LinearizedState(linear::MechanismState, current_state_vector::AbstractVector)
        nq = num_positions(linear)
        nv = num_velocities(linear)
        na = num_additional_states(linear)
        N = nq + nv + na
        xdiff = Vector{ForwardDiff.Dual{N, Float64}}(N)
        ForwardDiff.seed!(xdiff, state_vector(linear), ForwardDiff.construct_seeds(ForwardDiff.Partials{N, Float64}))
        diffstate = MechanismState(linear.mechanism, xdiff[1:nq], xdiff[nq + (1:nv)], xdiff[nq + nv + (1:na)])
        LinearizedState(linear, current_state_vector, diffstate)
    end

    unwrap(p::Point3D) = (v -> Point3D(p.frame, v), p.v)
    unwrap(p::FreeVector3D) = (v -> FreeVector3D(p.frame, v), p.v)
    unwrap(p) = (identity, p)

    configuration(s::LinearizedState, joint::Joint) = 
        @view s.current_state_vector[parentindexes(configuration(s.dual_state, joint))...]

    linearization_state(s::LinearizedState) = s.linearization_state
    linearization_state_vector(s::LinearizedState) = state_vector(s.linearization_state)

    function evaluate(f::Function, s::LinearizedState)
        wrapper, ydual = unwrap(f(s.dual_state))
        nx = length(s.current_state_vector)
        v = ForwardDiff.value.(ydual)
        Δx = s.current_state_vector .- linearization_state_vector(s)
        
        if isa(v, AbstractVector)
            J = similar(v, (length(v), nx))
            ForwardDiff.extract_jacobian!(J, ydual, nx)
            wrapper(v .+ J * Δx)
        else
            wrapper(v + ForwardDiff.partials(ydual)' * Δx)
        end 
    end

    function jacobian(f::Function, s::LinearizedState)
        wrapper, ydual = unwrap(f(s.dual_state))
        nx = length(s.current_state_vector)
        v = ForwardDiff.value.(ydual)
        J = similar(v, (length(v), nx))
        ForwardDiff.extract_jacobian!(J, ydual, nx)
        J
    end

end

using .Linear: LinearizedState, linearization_state

include("environments.jl")
include("conditional_jump_extensions.jl")
include("variable_containers.jl")
include("contact.jl")
include("joint_limits.jl")
include("simulation.jl")

end # module
