const MAX_JOINT_FORCE = 200

generalized_force(r::JointLimitResult) = r.λ * r.scaling

_lower(b::RigidBodyDynamics.Bounds) = b.lower
_upper(b::RigidBodyDynamics.Bounds) = b.upper

function resolve_joint_limit(model::Model, xnext::LinearizedState, joint::Joint)
    mechanism = xnext.linearization_state.mechanism
    total_weight = mass(mechanism) * norm(mechanism.gravitational_acceleration)
    bounds = RigidBodyDynamics.position_bounds(joint)
    N = length(bounds)
    if N > 0
        lb = _lower.(bounds)
        ub = _upper.(bounds)

        λ = @variable(model, [1:N],
                      lowerbound=-MAX_JOINT_FORCE,
                      upperbound=MAX_JOINT_FORCE,
                      basename="{λ_{joint}}")
        q = current_configuration(xnext, joint)

        @constraints model begin
            [i=1:N], q[i] >= lb[i]
            [i=1:N], q[i] <= ub[i]
        end

        for i in 1:length(λ)
            @disjunction(model, (q[i] >= ub[i]), (λ[i] >= 0))
            @disjunction(model, (q[i] <= lb[i]), (λ[i] <= 0))
        end

        JointLimitResult(λ, total_weight)
    else
        JointLimitResult(Vector{Variable}(), total_weight)
    end
end
