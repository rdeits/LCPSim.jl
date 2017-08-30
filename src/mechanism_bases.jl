function planar_base()
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world; gravity=SVector(0, 0, -9.81))

    frame = CartesianFrame3D("dummy")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    dummy = RigidBody(inertia)
    base_x = Joint("base_x", Prismatic([1., 0, 0]))
    position_bounds(base_x) .= Bounds(-10., 10)
    velocity_bounds(base_x) .= Bounds(-10., 10)
    effort_bounds(base_x) .= Bounds(0., 0)
    attach!(mechanism, world, base_x, eye(Transform3D, frame_before(base_x), default_frame(world)), dummy)

    frame = CartesianFrame3D("base")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    base = RigidBody(inertia)
    base_z = Joint("base_z", Prismatic([0., 0, 1]))
    position_bounds(base_z) .= Bounds(-10., 10)
    velocity_bounds(base_z) .= Bounds(-10., 10)
    effort_bounds(base_z) .= Bounds(0., 0)
    attach!(mechanism, dummy, base_z, eye(Transform3D, frame_before(base_z), default_frame(dummy)), base)
    mechanism, base
end

function planar_revolute_base()
    mechanism, base = planar_base()
    frame = CartesianFrame3D("base_revolute")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    joint = Joint("base_rotation", Revolute([0., 1., 0]))
    position_bounds(joint) .= Bounds(-π, π)
    velocity_bounds(joint) .= Bounds(-2π, 2π)
    effort_bounds(joint) .= Bounds(0., 0)
    attach!(mechanism, base, joint, eye(Transform3D, frame_before(joint), default_frame(base)), body)
    mechanism, body
end
