using JuMP: AbstractJuMPScalar

struct ContactResult{T, M}
    β::Vector{T}
    λ::T
    c_n::T
    Δr::FreeVector3D{SVector{3, T}}
    point::Point3D{SVector{3, M}}
    obs::Obstacle{M}
end

struct JointLimitResult{T}
    λ::Vector{T}
    # direction::Vector{M}
end

struct LCPUpdate{T, M, U}
    state::StateRecord{T, M}
    input::Vector{U}
    contacts::Dict{RigidBody{M}, Vector{ContactResult{T, M}}}
    joint_contacts::Vector{JointLimitResult{T}}
end

LCPUpdate(state::StateRecord{T}, input::AbstractVector{U}, contacts::Associative{<:RigidBody{M}}, joint_contacts::Associative) where {T, M, U} =
    LCPUpdate{T, M, U}(state, input, contacts, joint_contacts)

_getvalue(x::AbstractVector{<:Number}) = x
_getvalue(x::AbstractVector{<:AbstractJuMPScalar}) = getvalue(x)


JuMP.getvalue(r::JointLimitResult) = JointLimitResult(getvalue.(r.λ))
JuMP.getvalue(r::StateRecord{<:AbstractJuMPScalar}) = StateRecord(r.mechanism, getvalue.(r.state))
JuMP.getvalue(up::LCPUpdate) =
    LCPUpdate(getvalue(up.state),
              _getvalue(up.input),
              Dict([k => getvalue.(v) for (k, v) in up.contacts]),
              getvalue.(up.joint_contacts))
JuMP.getvalue(f::FreeVector3D) = FreeVector3D(f.frame, getvalue(f.v))
JuMP.getvalue(c::ContactResult) = ContactResult(getvalue.((c.β, c.λ, c.c_n, c.Δr))..., c.point, c.obs)

function JuMP.setvalue(r::JointLimitResult{<:AbstractJuMPScalar}, seed::JointLimitResult{<:Number})
    setvalue.(r.λ, seed.λ)
end

function JuMP.setvalue(up::LCPUpdate{<:AbstractJuMPScalar}, seed::LCPUpdate{<:Number})
    setvalue(up.state, seed.state)
    if eltype(up.input) <: AbstractJuMPScalar
        for i in eachindex(up.input)
            var = up.input[i]
            if var.m.colCat[var.col] != :Fixed
                setvalue(var, seed.input[i])
            end
        end
    end
    setvalue.(up.contacts, seed.contacts)
    setvalue.(up.joint_contacts, seed.joint_contacts)
end

JuMP.setvalue(f::FreeVector3D{<:AbstractVector{S}}, seed::FreeVector3D{<:AbstractVector{N}}) where {S <: AbstractJuMPScalar, N <: Number} =
    setvalue.(f.v, seed.v)

function JuMP.setvalue(contact::ContactResult{<:AbstractJuMPScalar}, seed::ContactResult{<:Number})
    @assert contact.obs == seed.obs
    setvalue(contact.β, seed.β)
    setvalue(contact.λ, seed.λ)
    setvalue(contact.c_n, seed.c_n)
    setvalue(contact.Δr, seed.Δr)
end

function JuMP.setvalue(contacts::Dict{<:RigidBody, <:Vector{<:ContactResult}}, seeds::Dict)
    for body in keys(contacts)
        setvalue.(contacts[body], seeds[body])
    end
end

function JuMP.setvalue(state::StateRecord{<:AbstractJuMPScalar}, seed::StateRecord{<:Number})
    setvalue.(state.state, seed.state)
end
