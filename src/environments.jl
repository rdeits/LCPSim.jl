
struct Obstacle{T}
    frame::CartesianFrame3D
    interior::SimpleHRepresentation{3, T}
    contact_face::HalfSpace{3, T}
    μ::T
end

struct FreeRegion{T}
    frame::CartesianFrame3D
    interior::SimpleHRepresentation{3, T}
end

struct ContactEnvironment{T}
    points::Vector{Point3D{SVector{3, T}}}
    obstacles::Vector{Obstacle{T}}
    free_regions::Vector{FreeRegion{T}}
end

struct Environment{T}
    contacts::Dict{RigidBody{T}, ContactEnvironment{T}}
end

function planar_obstacle(frame, normal::AbstractVector{T}, point::AbstractVector{T}, μ=1.0) where T
    normal = normalize(normal)
    b = normal' * point
    Obstacle(frame, 
             SimpleHRepresentation{3, T}(reshape(normal, (1, 3)), [b]),
             HalfSpace{3, T}(normal, b),
             μ)
end

function space_between(walls::Vector{<:Obstacle})
    frame = first(walls).frame
    @assert all([obs.frame == frame for obs in walls])
    @assert all([length(obs.interior) == 1 for obs in walls])

    FreeRegion(frame, SimpleHRepresentation(
        vcat([-obs.interior.A for obs in walls]...),
        vcat([-obs.interior.b for obs in walls]...)))
end

