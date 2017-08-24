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

include("linearized_mechanism.jl")
using .Linear: LinearizedState, linearization_state

include("environments.jl")
include("conditional_jump_extensions.jl")
include("variable_containers.jl")
include("contact.jl")
include("joint_limits.jl")
include("simulation.jl")

end # module
