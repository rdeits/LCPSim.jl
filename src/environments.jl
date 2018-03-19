
struct Obstacle{T}
    frame::CartesianFrame3D
    interior::SimpleHRepresentation{3, T}
    contact_face::HalfSpace{3, T}
    μ::T
    contact_basis::Vector{FreeVector3D{SVector{3, T}}}
end

function Obstacle(frame::CartesianFrame3D, interior::SimpleHRepresentation, contact_face::HalfSpace, μ, motion_type::Symbol)
    basis = contact_basis(frame, contact_face, μ, motion_type)
    Obstacle(frame, interior, contact_face, μ, basis)
end

contact_basis(obs::Obstacle) = obs.contact_basis

function contact_basis(frame::CartesianFrame3D, contact_face::HalfSpace{3, T}, μ, motion_type::Symbol) where T
    a::SVector{3,T} = normalize(SVector{3, T}(contact_face.a))
    if motion_type == :xz
        return [
            FreeVector3D(frame, Rotations.RotY(π/2) * a),
            FreeVector3D(frame, Rotations.RotY(-π/2) * a)
            ]
    elseif motion_type == :xy
        return [
            FreeVector3D(frame, Rotations.RotZ(π/2) * a),
            FreeVector3D(frame, Rotations.RotZ(-π/2) * a)
            ]
    elseif motion_type == :yz
        return [
            FreeVector3D(frame, Rotations.RotX(π/2) * a),
            FreeVector3D(frame, Rotations.RotX(-π/2) * a)
            ]
    elseif motion_type == :xyz
        R = Rotations.rotation_between(SVector{3, T}(0, 0, 1), a)
        D = [FreeVector3D(frame, R * RotZ(θ) * SVector(1, 0, 0))
                for θ in 0:π/2:3π/2]
        for i in 1:length(D)
            @assert isapprox(dot(D[i].v, contact_face.a), 0, atol=1e-15)
            if i > 1
                @assert isapprox(dot(D[i].v, D[i-1].v), 0, atol=1e-15)
            end
        end
        return D
    else
        throw(ArgumentError("Unrecognized motion type: $motion_type (should be :xy, :xz, :yz, or :xyz)"))
    end
end

struct Environment{T}
    contacts::Vector{Tuple{RigidBody{T}, Point3D{SVector{3, T}}, Obstacle{T}}}
end

function planar_obstacle(frame, normal::AbstractVector{T}, point::AbstractVector{T}, μ=1.0, motion_type::Symbol=:xyz) where T
    normal = normalize(normal)
    b = normal' * point
    Obstacle(frame, 
             SimpleHRepresentation{3, T}(reshape(normal, (1, 3)), [b]),
             HalfSpace{3, T}(normal, b),
             μ,
             motion_type)
end
