using JuMP: AbstractJuMPScalar

struct ContactResult{T, M}
    β::Vector{T}
    λ::T
    c_n::T
    # Δr::FreeVector3D{SVector{3, T}}
    point::Point3D{SVector{3, M}}
    obs::Obstacle{M}
    scaling::M
end

Base.show(io::IO, cr::ContactResult) = print(io, """ContactResult(β=$(cr.β), λ=$(cr.λ), c_n=$(cr.c_n), ...)""")


struct JointLimitResult{T, M}
    λ::Vector{T}
    scaling::M
end

struct LCPUpdate{T, M, U, Tdt}
    Δt::Tdt
    state::StateRecord{T, M}
    input::Vector{U}
    contacts::Dict{RigidBody{M}, Vector{ContactResult{T, M}}}
    joint_contacts::Vector{JointLimitResult{T, M}}
end

LCPUpdate(Δt::Tdt, ::StateRecord{T}, input::AbstractVector{U}, contacts::Associative{<:RigidBody{M}}, joint_contacts::Associative) where {T, M, U, Tdt} =
    LCPUpdate{T, M, U, Tdt}(Δt, state, input, contacts, joint_contacts)

_getvalue(x::Variable) = JuMP._getValue(x)
_getvalue(x::Number) = x
_getvalue(f::FreeVector3D) = FreeVector3D(f.frame, _getvalue.(f.v))
function _getvalue(a::AffExpr)
    ret = a.constant
    for it in 1:length(a.vars)
        ret += a.coeffs[it] * _getvalue(a.vars[it])
    end
    ret
end
function _getvalue(a::QuadExpr)
    ret = _getvalue(a.aff)
    for it in 1:length(a.qvars1)
        ret += a.qcoeffs[it] * _getvalue(a.qvars1[it]) * _getvalue(a.qvars2[it])
    end
    return ret
end


JuMP.getvalue(r::JointLimitResult) = JointLimitResult(_getvalue.(r.λ), r.scaling)
JuMP.getvalue(r::StateRecord{<:AbstractJuMPScalar}) = StateRecord(r.mechanism, _getvalue.(r.state))
JuMP.getvalue(up::LCPUpdate) =
    LCPUpdate(_getvalue(up.Δt),
              getvalue(up.state),
              _getvalue.(up.input),
              Dict([k => getvalue.(v) for (k, v) in up.contacts]),
              getvalue.(up.joint_contacts))
JuMP.getvalue(c::ContactResult) = ContactResult(_getvalue.(c.β),
                                                _getvalue(c.λ),
                                                _getvalue(c.c_n),
                                                # _getvalue(c.Δr),
                                                c.point,
                                                c.obs,
                                                c.scaling)

function JuMP.setvalue(r::JointLimitResult{<:AbstractJuMPScalar}, seed::JointLimitResult{<:Number})
    setvalue.(r.λ, seed.λ)
end

function JuMP.setvalue(up::LCPUpdate{<:AbstractJuMPScalar}, seed::LCPUpdate{<:Number})
    if typeof(up.Δt) <: AbstractJuMPScalar
        setvalue(up.Δt, seed.Δt)
    end
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
    # setvalue(contact.Δr, seed.Δr)
end

function JuMP.setvalue(contacts::Dict{<:RigidBody, <:Vector{<:ContactResult}}, seeds::Dict)
    for body in keys(contacts)
        setvalue.(contacts[body], seeds[body])
    end
end

function JuMP.setvalue(state::StateRecord{<:AbstractJuMPScalar}, seed::StateRecord{<:Number})
    setvalue.(state.state, seed.state)
end
