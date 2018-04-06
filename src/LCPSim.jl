__precompile__()

module LCPSim

export Obstacle,
       FreeRegion,
       ContactEnvironment,
       Environment,
       planar_obstacle,
       simulate,
       update,
       optimize,
       planar_base,
       planar_revolute_base,
       LinearizedState,
       linearized

using StaticArrays
using JuMP
using JuMP: GenericAffExpr
using ConditionalJuMP
using Base.Test
using RigidBodyDynamics
using RigidBodyDynamics: HalfSpace3D, separation
using Rotations
using CoordinateTransformations: transform_deriv
using ForwardDiff
using MechanismGeometries
using GeometryTypes: HyperSphere, origin

include("state_record.jl")
include("linearized_mechanism.jl")
using .Linear

include("environments.jl")
include("variable_containers.jl")
include("contact.jl")
include("joint_limits.jl")
include("simulation.jl")
include("mechanism_bases.jl")
include("contact_lqr.jl")

end # module
