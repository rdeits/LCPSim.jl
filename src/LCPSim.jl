__precompile__()

module LCPSim

export Obstacle,
       FreeRegion,
       ContactEnvironment,
       Environment,
       planar_obstacle,
       space_between,
       simulate,
       update,
       optimize,
       planar_base,
       planar_revolute_base,
       LinearizedState,
       linearized

using Polyhedra
using StaticArrays
using JuMP
using JuMP: GenericAffExpr
using ConditionalJuMP
using Base.Test
using RigidBodyDynamics
using Rotations
using ForwardDiff

include("state_record.jl")
include("linearized_mechanism.jl")
using .Linear

include("environments.jl")
include("conditional_jump_extensions.jl")
include("variable_containers.jl")
include("contact.jl")
include("joint_limits.jl")
include("simulation.jl")
include("mechanism_bases.jl")

end # module
