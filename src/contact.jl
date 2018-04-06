contact_normal(obs::Obstacle) = obs.contact_face.outward_normal

RigidBodyDynamics.separation(obs::Obstacle, p::Point3D) = separation(obs.contact_face, p)

function contact_force(r::ContactResult)
    n = contact_normal(r.obs)
    D = contact_basis(r.obs)
    @framecheck(n.frame, D[1].frame)
    v = r.c_n .* n.v
    for i in eachindex(r.β)
        v += r.β[i] .* D[i].v
    end
    FreeVector3D(n.frame, r.scaling * v)
end

function append(v::AbstractVector{T}, x::T) where T
    result = similar(v, length(v) + 1)
    for i in 1:length(v)
        result[i] = v[i]
    end
    result[end] = x
    result
end

function add_contact_constraints(model::Model, μ, β, λ, c_n, D_transpose_times_v, separation_from_obstacle)

    @constraints model begin
        separation_from_obstacle >= 0 # (7)
        μ * c_n - sum(β) >= 0 # (9)
    end

    for d in D_transpose_times_v
        @constraint model λ + d >= 0 # (8)
    end

    @disjunction(model, separation_from_obstacle == 0, c_n == 0)
    for j in 1:length(D_transpose_times_v)
        d = D_transpose_times_v[j]
        λ_plus_d = AffExpr(append(d.vars, λ), append(d.coeffs, 1.0), d.constant)
        @disjunction(model, (λ_plus_d == 0), (β[j] == 0)) # (11)
    end
    @disjunction(model, (μ * c_n - sum(β) == 0), (λ == 0)) # (12)
end

function resolve_contact(xnext::LinearizedState, body::RigidBody, point::Point3D, obstacle::Obstacle, model::Model)
    mechanism = xnext.linearization_state.mechanism
    total_weight = mass(mechanism) * norm(mechanism.gravitational_acceleration)

    D = contact_basis(obstacle)
    β   = @variable(model, [1:length(D)], lowerbound=0, basename="β", upperbound=100)
    λ   = @variable(model, lowerbound=0, basename="λ", upperbound=100)
    c_n = @variable(model, lowerbound=0, basename="c_n", upperbound=100)
    world = default_frame(root_body(linearization_state(xnext).mechanism))

    separation_from_obstacle = linearized(xnext) do x
        separation(obstacle, transform(x, point, obstacle.contact_face.outward_normal.frame))
    end


    D_transpose_times_v = map(D) do d
        linearized(xnext) do x
            v = point_velocity(twist_wrt_world(x, body), transform_to_root(linearization_state(xnext), point.frame) * point)
            dot(transform_to_root(linearization_state(xnext), d.frame) * d,
                v)
        end
    end

    add_contact_constraints(model, obstacle.μ, β, λ, c_n, D_transpose_times_v, separation_from_obstacle)
    ContactResult(β, λ, c_n, point, obstacle, total_weight)
end

# function add_contact_constraints_nonsliding(model::Model, point, obstacle, β, λ, c_n, D, separation_from_obstacle, contact_velocity)
#     D_transpose_times_v = [dot(d, contact_velocity) for d in D]

#     @constraints model begin
#         separation_from_obstacle >= 0 # (7)
#         # λ .+ D_transpose_times_v .>= 0 # (8)
#         obstacle.μ * c_n .- sum(β) >= 0 # (9)
#     end

#     @disjunction(model,
#                  (&)(separation_from_obstacle == 0,
#                      [@?(contact_velocity.v[i] == 0) for i in 1:3]...),
#                  c_n == 0)
#     ContactResult(β, λ, c_n, point, obstacle)
# end

# function add_contact_constraints_sticking(model::Model, point, obstacle, β, λ, c_n, D, separation_from_obstacle, contact_velocity)

#     @constraints model begin
#         separation_from_obstacle <= 1e-3 # (7)
#         separation_from_obstacle >= -1e-3
#         contact_velocity.v .== 0
#         λ == 0 # (8)
#         obstacle.μ * c_n .- sum(β) >= 0 # (9)
#     end

#     ContactResult(β, λ, c_n, point, obstacle)
# end

