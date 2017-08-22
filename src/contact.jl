
function contact_basis(obs::Obstacle)
    θ = atan(obs.μ)
    R = RotY(θ)
    SVector(
        FreeVector3D(obs.frame, R * obs.contact_face.a), 
        FreeVector3D(obs.frame, R' * obs.contact_face.a))
end

contact_normal(obs::Obstacle) = FreeVector3D(obs.frame, obs.contact_face.a)

function separation(obs::Obstacle, p::Point3D)
    @framecheck obs.frame p.frame
    n = contact_normal(obs)
    n.v' * p.v - obs.contact_face.β
end

_vec(f::FreeVector3D) = convert(Vector, f.v)
_vec(p::Point3D) = convert(Vector, p.v)

function contact_force(r::ContactResult)
    n = contact_normal(r.obs)
    D = contact_basis(r.obs)
    @framecheck(n.frame, D[1].frame)
    # r.c_n * n .+ D * r.β
    FreeVector3D(n.frame, r.c_n .* n.v .+ sum(broadcast.(*, _vec.(D), r.β)))
end


function resolve_contact(xnext::MechanismState, body::RigidBody, point::Point3D, obstacle::Obstacle, model::Model, x_dynamics::MechanismState{<:Number})
    D = contact_basis(obstacle)
    k = length(D)

    β = @variable(model,   [1:k], lowerbound=0,   basename="β",     upperbound=100)
    λ = @variable(model,          lowerbound=0,   basename="λ",     upperbound=100)
    c_n = @variable(model,        lowerbound=0,   basename="c_n",   upperbound=100)

    function _separation(x)
        q = x[1:num_positions(xnext)]
        v = x[(num_positions(xnext)+1):end]
        x_diff = MechanismState(xnext.mechanism, q, v)
        point_in_world = transform_to_root(x_diff, point.frame) * point
        separation(obstacle, point_in_world)
    end

    separation_from_obstacle = _separation(state_vector(x_dynamics)) + (state_vector(xnext) - state_vector(x_dynamics))' * ForwardDiff.gradient(_separation, state_vector(x_dynamics))

    function _contact_velocity(x)
        q = x[1:num_positions(xnext)]
        v = x[(num_positions(xnext)+1):end]
        x_diff = MechanismState(xnext.mechanism, q, v)
        point_in_world = transform_to_root(x_diff, point.frame) * point
        contact_velocity = point_velocity(twist_wrt_world(x_diff, body), point_in_world).v
    end

    contact_velocity = FreeVector3D(root_frame(xnext.mechanism), _contact_velocity(state_vector(x_dynamics)) + ForwardDiff.jacobian(_contact_velocity, state_vector(x_dynamics)) * (state_vector(xnext) - state_vector(x_dynamics)))

    D_transpose_times_v = [dot(d, contact_velocity) for d in D]

    @constraints model begin
        λ .+ D_transpose_times_v .>= 0 # (8)
        obstacle.μ * c_n .- sum(β) >= 0 # (9)
    end

    @disjunction(model, (separation_from_obstacle == 0), (c_n == 0)) # (10)
    for j in 1:k
        @disjunction(model, ((λ + D_transpose_times_v[j]) == 0), β[j] == 0) # (11)
    end
    @disjunction(model, (obstacle.μ * c_n - sum(β) == 0), (λ == 0)) # (12)

    ContactResult(β, λ, c_n, point, obstacle)
end
