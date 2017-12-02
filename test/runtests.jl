using LCPSim
using Base.Test

# write your own tests here
include("contact_lqr.jl")
include("linearized_mechanism.jl")
include("warmstart.jl")
include("inclined_brick.jl")
include("planar_vs_dummy_link.jl")
include("acrobot.jl")
if Pkg.installed("Gurobi") != nothing && Pkg.installed("RigidBodyTreeInspector") != nothing
    include("../examples/box.jl")
end
