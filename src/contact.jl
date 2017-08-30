
function contact_basis(obs::Obstacle{T}) where T
    θ = π/2
    R = RotY(θ)
    a::SVector{3,T} = SVector{3, T}(obs.contact_face.a)
    SVector(
        FreeVector3D(obs.frame, R * a),
        FreeVector3D(obs.frame, R' * a))
end

contact_normal(obs::Obstacle{T}) where {T} = 
    FreeVector3D(obs.frame, (SVector{3, T}(obs.contact_face.a)::SVector{3, T}))

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

function add_contact_constraints(model::Model, point, obstacle, β, λ, c_n, D, separation_from_obstacle, contact_velocity)
    D_transpose_times_v = [dot(d, contact_velocity) for d in D]
    k = length(D)

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

function resolve_contact(xnext::LinearizedState, body::RigidBody, point::Point3D, obstacle::Obstacle, model::Model)
    D = contact_basis(obstacle)
    k = length(D)

    β   = @variable(model,   [1:k], lowerbound=0,   basename="β",     upperbound=200)
    λ   = @variable(model,          lowerbound=0,   basename="λ",     upperbound=200)
    c_n = @variable(model,        lowerbound=0,   basename="c_n",   upperbound=200)

    point_in_world = linearized(x -> transform_to_root(x, point.frame) * point, xnext)
    separation_from_obstacle = separation(obstacle, point_in_world)

    contact_velocity = linearized(x -> point_velocity(twist_wrt_world(x, body), transform_to_root(linearization_state(xnext), point.frame) * point), xnext)
    add_contact_constraints(model, point, obstacle, β, λ, c_n, D, separation_from_obstacle, contact_velocity)
end
