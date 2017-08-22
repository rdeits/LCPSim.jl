
struct ContactResult{T, M}
    β::Vector{T}
    λ::T
    c_n::T
    point::Point3D{SVector{3, M}}
    obs::Obstacle{M}
end

struct JointLimitResult{T, M}
    λ::T
    direction::Vector{M}
end

struct LCPUpdate{T, M, S <: MechanismState{T, M}, I <:AbstractVector}
    state::S
    input::I
    contacts::Dict{RigidBody{M}, Vector{ContactResult{T, M}}}
    joint_contacts::Dict{Joint, Vector{JointLimitResult{T, M}}}
end

LCPUpdate(state::S, input::I, contacts::Associative, joint_contacts::Associative) where {T, M, S <: MechanismState{T, M}, I <: AbstractVector} = 
    LCPUpdate{T, M, S, I}(state, input, contacts, joint_contacts)

_getvalue(x::AbstractVector{<:Number}) = x
_getvalue(x::AbstractVector{<:JuMP.AbstractJuMPScalar}) = getvalue(x)


JuMP.getvalue(r::JointLimitResult) = JointLimitResult(getvalue(r.λ), r.direction)
JuMP.getvalue(x::MechanismState{<:JuMP.AbstractJuMPScalar}) = 
    MechanismState(x.mechanism, getvalue(configuration(x)), getvalue(velocity(x)), getvalue(additional_state(x)))
# JuMP.getvalue(p::Pair{<:RigidBody, <:AbstractVector{<:ContactResult}}) = p.first => getvalue.(p.second)
# JuMP.getvalue(p::Pair{<:Joint, <:AbstractVector{<:JointLimitResult}}) = p.first => getvalue.(p.second)
JuMP.getvalue(up::LCPUpdate) =
    LCPUpdate(getvalue(up.state), 
              _getvalue(up.input), 
              Dict([k => getvalue.(v) for (k, v) in up.contacts]),
              Dict([k => getvalue.(v) for (k, v) in up.joint_contacts]))
              # map(getvalue, up.contacts), map(getvalue, up.joint_contacts))
JuMP.getvalue(c::ContactResult) = ContactResult(getvalue.((c.β, c.λ, c.c_n))..., c.point, c.obs)
JuMP.getvalue(f::FreeVector3D) = FreeVector3D(f.frame, getvalue(f.v))

function JuMP.setvalue(r::JointLimitResult{<:JuMP.AbstractJuMPScalar}, seed::JointLimitResult{<:Number})
    setvalue(r.λ, seed.λ)
    @assert getvalue(r.generalized_force) ≈ seed.generalized_force
end

function JuMP.setvalue(d::Dict{<:Joint, <:Vector{<:JointLimitResult}}, seed::Dict)
    [setvalue.(v1, v2) for (v1, v2) in zip(values(d), values(seed))]
end

function JuMP.setvalue(up::LCPUpdate{<:JuMP.AbstractJuMPScalar}, seed::LCPUpdate{<:Number})
    setvalue(up.state, seed.state)
    setvalue.(up.contacts, seed.contacts)
    setvalue(up.joint_contacts, seed.joint_contacts)
end

function JuMP.setvalue(contact::ContactResult{<:JuMP.AbstractJuMPScalar}, seed::ContactResult{<:Number})
    @assert contact.obs == seed.obs
    setvalue(contact.β, seed.β)
    setvalue(contact.λ, seed.λ)
    setvalue(contact.c_n, seed.c_n)
end
