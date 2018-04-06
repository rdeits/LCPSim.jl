
struct Obstacle{T}
    interior::Vector{HalfSpace3D{T}}
    contact_face::HalfSpace3D{T}
    μ::T
    contact_basis::Vector{FreeVector3D{SVector{3, T}}}
end

function Obstacle(interior::AbstractVector{<:HalfSpace3D}, contact_face::HalfSpace3D, μ, motion_type::Symbol)
    basis = contact_basis(contact_face, μ, motion_type)
    Obstacle(interior, contact_face, μ, basis)
end

contact_basis(obs::Obstacle) = obs.contact_basis

function contact_basis(contact_face::HalfSpace3D{T}, μ, motion_type::Symbol) where T
    a = contact_face.outward_normal.v
    frame = contact_face.outward_normal.frame
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
            @assert isapprox(dot(D[i].v, a), 0, atol=1e-15)
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

function planar_obstacle(normal::FreeVector3D, point::Point3D, μ=1.0, motion_type::Symbol=:xyz)
    normal = normalize(normal)
    face = HalfSpace3D(point, normal)
    Obstacle([face],
             face,
             μ,
             motion_type)
end

planar_obstacle(frame::CartesianFrame3D, normal::AbstractVector, point::AbstractVector, args...) =
    planar_obstacle(FreeVector3D(frame, normal), Point3D(frame, point), args...)