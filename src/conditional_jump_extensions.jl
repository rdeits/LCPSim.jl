Base.@pure ConditionalJuMP.isjump(::Point3D{<:AbstractArray{<:JuMP.AbstractJuMPScalar}}) = true

function ConditionalJuMP.Conditional(op::typeof(in), x::Point3D, P::FreeRegion)
    @framecheck(x.frame, P.frame)
    ConditionalJuMP.Conditional(&, [@?(P.interior.A[i, :]' * x.v <= P.interior.b[i]) for i in 1:length(P.interior)]...)
end
