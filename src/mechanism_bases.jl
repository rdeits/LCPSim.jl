function planar_base()
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world; gravity=SVector(0, 0, -9.81))

    frame = CartesianFrame3D("dummy")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    dummy = RigidBody(inertia)
    base_x = Joint("base_x", Prismatic([1., 0, 0]))
    position_bounds(base_x) .= Bounds(-10., 10)
    velocity_bounds(base_x) .= Bounds(-1000., 1000)
    effort_bounds(base_x) .= Bounds(0., 0)
    attach!(mechanism, world, dummy, base_x)

    frame = CartesianFrame3D("base")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    base = RigidBody(inertia)
    base_z = Joint("base_z", Prismatic([0., 0, 1]))
    position_bounds(base_z) .= Bounds(-10., 10)
    velocity_bounds(base_z) .= Bounds(-1000., 1000)
    effort_bounds(base_z) .= Bounds(0., 0)
    attach!(mechanism, dummy, base, base_z)
    mechanism, base
end

function planar_revolute_base()
    mechanism, base = planar_base()
    frame = CartesianFrame3D("base_revolute")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    joint = Joint("base_rotation", Revolute([0., 1., 0]))
    position_bounds(joint) .= Bounds(-2π, 2π)
    velocity_bounds(joint) .= Bounds(-1000., 1000)
    effort_bounds(joint) .= Bounds(0., 0)
    attach!(mechanism, base, body, joint)
    mechanism, body
end

function quaternion_floating_base()
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world; gravity=SVector(0, 0, -9.81))
    frame = CartesianFrame3D("base")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    joint = Joint("floating_base", QuaternionFloating{Float64}())
    position_bounds(joint) .= Bounds(-10, 10)
    velocity_bounds(joint) .= Bounds(-1000, 1000)
    attach!(mechanism, world, body, joint)
    mechanism, body
end

function xyz_rpy_floating_base()
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world; gravity=SVector(0, 0, -9.81))

    parent = world
    frame = CartesianFrame3D("dummy_x")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    # joint = Joint("base_x", Prismatic([1., 0, 0]))
    joint = Joint("base_x", Fixed{Float64}())
    attach!(mechanism, parent, body, joint)

    parent = body
    frame = CartesianFrame3D("dummy_y")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    joint = Joint("base_y", Prismatic([0., 1, 0]))
    attach!(mechanism, parent, body, joint)

    parent = body
    frame = CartesianFrame3D("dummy_z")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    joint = Joint("base_z", Prismatic([0., 0, 1]))
    attach!(mechanism, parent, body, joint)

    parent = body
    frame = CartesianFrame3D("dummy_rx")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    joint = Joint("base_rx", Revolute([1., 0, 0]))
    attach!(mechanism, parent, body, joint)

    parent = body
    frame = CartesianFrame3D("dummy_ry")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    # joint = Joint("base_ry", Revolute([0., 1, 0]))
    joint = Joint("base_ry", Fixed{Float64}())
    attach!(mechanism, parent, body, joint)

    parent = body
    frame = CartesianFrame3D("dummy_rz")
    inertia = SpatialInertia(frame, 0 * eye(3), zeros(3), 0.0)
    body = RigidBody(inertia)
    joint = Joint("base_rz", Revolute([0., 0, 1]))
    attach!(mechanism, parent, body, joint)

    for joint in joints(mechanism)
        position_bounds(joint) .= Bounds(-10, 10)
        velocity_bounds(joint) .= Bounds(-1000, 1000)
        effort_bounds(joint) .= Bounds(0, 0)
    end

    mechanism, body
end
