
function resolve_joint_limit(xnext::MechanismState, joint::Joint, a::AbstractVector, b::Number, model::Model)
    λ = @variable(model, lowerbound=0, upperbound=100, basename="λ")
    q = configuration(xnext, joint)
    separation = a' * q - b
    @constraint model separation <= 0
    @disjunction(model, separation == 0, λ == 0)

    JointLimitResult(λ, -a)
end

function resolve_joint_limits(xnext::MechanismState, joint::Joint, limits::HRepresentation, model::Model)
    [resolve_joint_limit(xnext, joint, limits.A[i, :], limits.b[i], model) for i in 1:length(limits)]
end

generalized_force(r::JointLimitResult) = r.λ * r.direction
