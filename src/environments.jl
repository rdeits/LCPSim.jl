
struct Obstacle{T}
    frame::CartesianFrame3D
    interior::SimpleHRepresentation{3, T}
    contact_face::HalfSpace{3, T}
    μ::T
end

struct ContactEnvironment{T}
    points::Vector{Point3D{SVector{3, T}}}
    obstacles::Vector{Obstacle{T}}
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
