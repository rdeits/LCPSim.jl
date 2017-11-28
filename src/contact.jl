
function contact_basis(obs::Obstacle{T}) where T
    a::SVector{3,T} = normalize(SVector{3, T}(obs.contact_face.a))
    R = Rotations.rotation_between(SVector{3, T}(0, 0, 1), a)
    D = SVector{4}([
        FreeVector3D(obs.frame, R * RotZ(θ) * SVector(1, 0, 0))
        for θ in 0:π/2:3π/2])
    for i in 1:length(D)
        @assert isapprox(dot(D[i].v, obs.contact_face.a), 0, atol=1e-15)
        if i > 1
            @assert isapprox(dot(D[i].v, D[i-1].v), 0, atol=1e-15)
        end
    end
    D
end

contact_normal(obs::Obstacle{T}) where {T} = 
    FreeVector3D(obs.frame, (SVector{3, T}(obs.contact_face.a)::SVector{3, T}))

function separation(obs::Obstacle, p::Point3D)
    @framecheck obs.frame p.frame
    n = contact_normal(obs)
    n.v' * p.v - obs.contact_face.β
end

function contact_force(r::ContactResult{<:JuMP.AbstractJuMPScalar})
    n = contact_normal(r.obs)
    D = contact_basis(r.obs)
    @framecheck(n.frame, D[1].frame)
    v = r.c_n .* n.v
    for i in eachindex(r.β)
        di = D[i].v
        for j in 1:length(di)
            append!(v[j], r.β[i] * di[j])
        end
    end
    FreeVector3D(n.frame, v)
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

function append(v::AbstractVector{T}, x::T) where T
    result = similar(v, length(v) + 1)
    for i in 1:length(v)
        result[i] = v[i]
    end
    result[end] = x
    result
end

function add_contact_constraints(model::Model, point, obstacle, β, λ, c_n, D_transpose_times_v, separation_from_obstacle, Δr, Δr_true)

    @constraints model begin
        separation_from_obstacle >= 0 # (7)
        obstacle.μ * c_n - sum(β) >= 0 # (9)
    end

    for d in D_transpose_times_v
        @constraint model λ + d >= 0 # (8)
    end

    @disjunction(model, 
        ((separation_from_obstacle <= 1e-3) & (Δr.v .== Δr_true.v)),
        ((c_n == 0) & (Δr.v .== 0))
    ) # (10)
    for j in 1:length(D_transpose_times_v)
        d = D_transpose_times_v[j]
        λ_plus_d = AffExpr(append(d.vars, λ), append(d.coeffs, 1.0), d.constant)
        @disjunction(model, (λ_plus_d == 0), (β[j] == 0)) # (11)
    end
    @disjunction(model, (obstacle.μ * c_n - sum(β) == 0), (λ == 0)) # (12)

    ContactResult(β, λ, c_n, Δr, point, obstacle)
end

function resolve_contact(xnext::LinearizedState, body::RigidBody, point::Point3D, obstacle::Obstacle, model::Model)
    D = contact_basis(obstacle)
    β   = @variable(model, [1:length(D)], lowerbound=0, basename="β", upperbound=1000)
    λ   = @variable(model, lowerbound=0, basename="λ", upperbound=1000)
    c_n = @variable(model, lowerbound=0, basename="c_n", upperbound=1000)
    world = default_frame(root_body(linearization_state(xnext).mechanism))
    Δr = FreeVector3D(world, SVector{3, Variable}(@variable(model, [1:3], lowerbound=-10, upperbound=10, basename="Δr")))

    Δr_true = linearized(xnext) do x
        (transform_to_root(x, point.frame) * point - center_of_mass(x))
    end - (transform_to_root(linearization_state(xnext), point.frame) * point - center_of_mass(linearization_state(xnext)))

    separation_from_obstacle = linearized(xnext) do x
        separation(obstacle, transform_to_root(x, point.frame) * point)
    end

    # contact_velocity = linearized(x -> point_velocity(twist_wrt_world(x, body), transform_to_root(linearization_state(xnext), point.frame) * point), xnext)
    # D_transpose_times_v = map(D) do d
    #     dot(d, contact_velocity)
    # end

    D_transpose_times_v = [
        linearized(x -> dot(d, point_velocity(twist_wrt_world(x, body), transform_to_root(linearization_state(xnext), point.frame) * point)), xnext) for d in D
        ]

    add_contact_constraints(model, point, obstacle, β, λ, c_n, D_transpose_times_v, separation_from_obstacle, Δr, Δr_true)
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

