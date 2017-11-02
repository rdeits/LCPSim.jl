
const VectorView{T} = SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
@assert typeof(view(zeros(5), 1:3)) === VectorView{Float64}

struct StateRecord{T, M}
    mechanism::Mechanism{M}
    state::Vector{T}
    configuration::VectorView{T}
    velocity::VectorView{T}
    additional_state::VectorView{T}

    function StateRecord{T, M}(mechanism::Mechanism{M}, state::AbstractVector{T}) where {T, M}
        c = view(state, 1:num_positions(mechanism))
        v = view(state, num_positions(mechanism) + (1:num_velocities(mechanism)))
        a = view(state, num_positions(mechanism) + num_velocities(mechanism) + (1:num_additional_states(mechanism)))
        new{T, M}(mechanism, state, c, v, a)
    end

    StateRecord{T}(mechanism::Mechanism{M}, state::AbstractVector{T}) where {T, M} = 
        StateRecord{T, M}(mechanism, state)
end

StateRecord(mechanism::Mechanism{M}, state::AbstractVector{T}) where {T, M} =
    StateRecord{T, M}(mechanism, state)

Base.convert(::Type{<:StateRecord{T}}, x::MechanismState) where {T} = StateRecord{T}(x.mechanism, copy(state_vector(x)))
Base.convert(::Type{StateRecord}, x::MechanismState{T}) where {T} = convert(StateRecord{T}, x)
RigidBodyDynamics.state_vector(r::StateRecord) = r.state
RigidBodyDynamics.configuration(r::StateRecord) = r.configuration
RigidBodyDynamics.velocity(r::StateRecord) = r.velocity
RigidBodyDynamics.set_configuration!(r::StateRecord, q::AbstractVector) = r.configuration .= q
RigidBodyDynamics.set_velocity!(r::StateRecord, v::AbstractVector) = r.velocity .= v
RigidBodyDynamics.num_positions(r::StateRecord) = length(r.configuration)
RigidBodyDynamics.num_velocities(r::StateRecord) = length(r.velocity)
