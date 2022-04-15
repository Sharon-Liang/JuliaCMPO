#module rrule
"""
    rrule for logtrexp function, eltype(M) <: Real
"""
function ChainRules.rrule(::typeof(logtrexp), M::AbstractArray{T}) where T<:Real
    e, v = symeigen(M)
    y = logsumexp(e)
    function logtrexp_pullback(ȳ)
        ∂y_∂M = v * diagm(exp.(e .- y)) * v' |> transpose
        M̄ = ȳ * ∂y_∂M
        return ChainRules.NoTangent(), M̄
    end
    return y, logtrexp_pullback
end

#end #module rrule

