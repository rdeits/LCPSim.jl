
function contact_basis(obs::Obstacle{T}) where T
    R = RotY(π/2)
    a::SVector{3,T} = SVector{3, T}(obs.contact_face.a)
    SVector{4}([
        FreeVector3D(obs.frame, RotZ(θ) * R * a)
        for θ in 0:π/2:3π/2])
end

contact_normal(obs::Obstacle{T}) where {T} = 
    FreeVector3D(obs.frame, (SVector{3, T}(obs.contact_face.a)::SVector{3, T}))

function separation(obs::Obstacle, p::Point3D)
    @framecheck obs.frame p.frame
    n = contact_normal(obs)
    n.v' * p.v - obs.contact_face.β
end

function contact_force(r::ContactResult)
    n = contact_normal(r.obs)
    D = contact_basis(r.obs)
    @framecheck(n.frame, D[1].frame)
    v = r.c_n .* n.v
    for i in eachindex(r.β)
        v += r.β[i] .* D[i].v
    end
    FreeVector3D(n.frame, v)
end

function add_contact_constraints(model::Model, point, obstacle, β, λ, c_n, D, separation_from_obstacle, contact_velocity)
    D_transpose_times_v = [dot(d, contact_velocity) for d in D]

    @constraints model begin
        separation_from_obstacle >= 0 # (7)
        λ .+ D_transpose_times_v .>= 0 # (8)
        obstacle.μ * c_n .- sum(β) >= 0 # (9)
    end

    # k = length(D)
    # cases = [
    #     @?(c_n == 0);
    #     [(&)(@?(separation_from_obstacle == 0), @?(λ + D_transpose_times_v[j] == 0), [@?(β[i] == 0) for i in 1:k if i != j]...) for j in 1:k]]
    # ConditionalJuMP.disjunction!(model, cases)

    @disjunction(model, (separation_from_obstacle <= 1e-3), (c_n == 0)) # (10)
    for j in 1:length(D)
        @disjunction(model, ((λ + D_transpose_times_v[j]) == 0), β[j] == 0) # (11)
    end
    @disjunction(model, (obstacle.μ * c_n - sum(β) == 0), (λ == 0)) # (12)

    ContactResult(β, λ, c_n, point, obstacle)
end

function resolve_contact(xnext::LinearizedState, body::RigidBody, point::Point3D, obstacle::Obstacle, model::Model)
    D = contact_basis(obstacle)
    β   = @variable(model, [1:length(D)], lowerbound=0, basename="β", upperbound=1000)
    λ   = @variable(model, lowerbound=0, basename="λ", upperbound=1000)
    c_n = @variable(model, lowerbound=0, basename="c_n", upperbound=1000)

    point_in_world = linearized(x -> transform_to_root(x, point.frame) * point, xnext)
    separation_from_obstacle = separation(obstacle, point_in_world)

    contact_velocity = linearized(x -> point_velocity(twist_wrt_world(x, body), transform_to_root(linearization_state(xnext), point.frame) * point), xnext)
    add_contact_constraints(model, point, obstacle, β, λ, c_n, D, separation_from_obstacle, contact_velocity)
end

function add_contact_constraints_nonsliding(model::Model, point, obstacle, β, λ, c_n, D, separation_from_obstacle, contact_velocity)
    D_transpose_times_v = [dot(d, contact_velocity) for d in D]

    @constraints model begin
        separation_from_obstacle >= 0 # (7)
        # λ .+ D_transpose_times_v .>= 0 # (8)
        obstacle.μ * c_n .- sum(β) >= 0 # (9)
    end

    @disjunction(model, 
                 (&)(separation_from_obstacle == 0, 
                     [@?(contact_velocity.v[i] == 0) for i in 1:3]...),
                 c_n == 0)
    ContactResult(β, λ, c_n, point, obstacle)
end

function add_contact_constraints_sticking(model::Model, point, obstacle, β, λ, c_n, D, separation_from_obstacle, contact_velocity)

    @constraints model begin
        separation_from_obstacle <= 1e-3 # (7)
        separation_from_obstacle >= -1e-3
        contact_velocity.v .== 0
        λ == 0 # (8)
        obstacle.μ * c_n .- sum(β) >= 0 # (9)
    end

    ContactResult(β, λ, c_n, point, obstacle)
end

