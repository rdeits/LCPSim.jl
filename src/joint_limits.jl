
function resolve_joint_limit(xnext::LinearizedState, joint::Joint, a::AbstractVector, b::Number, model::Model)
    λ = @variable(model, lowerbound=0, upperbound=200, basename="λ")
    q = current_configuration(xnext, joint)
    separation = a' * q - b
    @constraint model separation <= 0
    @disjunction(model, separation == 0, λ == 0)

    JointLimitResult(λ, -a)
end

function resolve_joint_limits(xnext::LinearizedState, joint::Joint, model::Model)
    A = vcat(eye(num_positions(joint)), -eye(num_positions(joint)))
    b = vcat(upper.(position_bounds(joint)), .-lower.(position_bounds(joint)))
    [resolve_joint_limit(xnext, joint, A[i, :], b[i], model) for i in 1:size(A, 1)]
end

generalized_force(r::JointLimitResult) = r.λ * r.direction
