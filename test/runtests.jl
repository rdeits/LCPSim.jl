using LCPSim
using Base.Test

@testset "LCPSim" begin
    include("planar_vs_dummy_link.jl")
    include("joint_limits.jl")
    include("contact_lqr.jl")
    include("linearized_mechanism.jl")
    include("warmstart.jl")
    include("inclined_brick.jl")
    include("acrobot.jl")
    include("../examples/box.jl")
    include("../examples/rimless_wheel.jl")
end
