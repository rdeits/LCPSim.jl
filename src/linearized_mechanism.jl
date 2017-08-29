
module Linear

    using RigidBodyDynamics
    import RigidBodyDynamics: configuration, velocity
    using ForwardDiff
    using MechanismComplementarity: StateRecord

    export set_current_configuration!,
           set_current_velocity!,
           set_linearization_configuration!,
           set_linearization_velocity!,
           current_configuration,
           current_velocity,
           LinearizedState,
           linearization_state

    struct LinearizedState{T, SCurrent <: StateRecord{T}, SLinear <: MechanismState{<:Number}, SDual <: MechanismState{<:ForwardDiff.Dual}}
        current_state::SCurrent
        linearization_state::SLinear
        dual_state::SDual

        function LinearizedState{T}(linear::S) where {T, S <: MechanismState{<:Number}}
            mechanism = linear.mechanism
            nq = num_positions(mechanism)
            nv = num_velocities(mechanism)
            na = num_additional_states(mechanism)
            N = nq + nv + na
            xdiff = Vector{ForwardDiff.Dual{N, Float64}}(N)
            ForwardDiff.seed!(xdiff, state_vector(linear), ForwardDiff.construct_seeds(ForwardDiff.Partials{N, Float64}))
            diffstate = MechanismState(mechanism, xdiff[1:nq], xdiff[nq + (1:nv)], xdiff[nq + nv + (1:na)])
            current = StateRecord(mechanism, Vector{T}(N))
            new{T, typeof(current), S, typeof(diffstate)}(StateRecord(mechanism, Vector{T}(N)),
                                         linear,
                                         diffstate)
        end
    end

    function LinearizedState(linear::MechanismState, current_state::StateRecord{T}) where T
        state = LinearizedState{T}(linear)
        set_current_configuration!(state, configuration(current_state))
        set_current_velocity!(state, velocity(current_state))
        state
    end

    set_current_configuration!(state::LinearizedState, q::AbstractVector) = set_configuration!(state.current_state, q)
    set_current_velocity!(state::LinearizedState, v::AbstractVector) = set_velocity!(state.current_state, v)
    function set_linearization_configuration!(state::LinearizedState, q::AbstractVector)
        set_configuration!(state.linearization_state, q)
        qdual = configuration(state.dual_state)
        qdual .= ForwardDiff.Dual.(q, ForwardDiff.partials.(qdual))
        setdirty!(state.dual_state)
    end
    function set_linearization_velocity!(state::LinearizedState, v::AbstractVector)
        set_velocity!(state.linearization_state, v)
        vdual = velocity(state.dual_state)
        vdual .= ForwardDiff.Dual.(v, ForwardDiff.partials.(vdual))
        setdirty!(state.dual_state)
    end

    unwrap(p::Point3D) = (v -> Point3D(p.frame, v), p.v)
    unwrap(p::FreeVector3D) = (v -> FreeVector3D(p.frame, v), p.v)
    unwrap(p::Transform3D) = (v -> Transform3D(p.from, p.to, v), p.mat)
    unwrap(p) = (identity, p)

    current_configuration(s::LinearizedState, joint::Joint) = 
        @view configuration(current_state(s))[parentindexes(configuration(s.dual_state, joint))...]
    current_velocity(s::LinearizedState, joint::Joint) = 
        @view velocity(current_state(s))[parentindexes(velocity(s.dual_state, joint))...]


    linearization_state(s::LinearizedState) = s.linearization_state
    linearization_state_vector(s::LinearizedState) = state_vector(s.linearization_state)
    current_state(s::LinearizedState) = s.current_state

    function evaluate(f::Function, s::LinearizedState)
        wrapper, ydual = unwrap(f(s.dual_state))
        nx = length(state_vector(current_state(s)))
        v = ForwardDiff.value.(ydual)
        Δx = state_vector(current_state(s)) .- state_vector(linearization_state(s))
        
        if isa(v, AbstractArray)
            J = similar(v, (length(v), nx))
            ForwardDiff.extract_jacobian!(J, ydual, nx)
            wrapper(v .+ reshape(J * Δx, size(v)))
        else
            wrapper(v + ForwardDiff.partials(ydual)' * Δx)
        end 
    end

    function jacobian(f::Function, s::LinearizedState)
        wrapper, ydual = unwrap(f(s.dual_state))
        nx = length(state_vector(current_state(s)))
        v = ForwardDiff.value.(ydual)
        J = similar(v, (length(v), nx))
        ForwardDiff.extract_jacobian!(J, ydual, nx)
        J
    end
end
