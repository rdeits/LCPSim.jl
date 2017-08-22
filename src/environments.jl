
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
