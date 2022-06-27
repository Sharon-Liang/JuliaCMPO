#module rrule
Zygote.@adjoint Array(x::CuArray) = Array(x), dy->(CuArray(dy),)

"""
    rrule for logtrexp function, eltype(M) <: Real
"""    
function ChainRules.rrule(::typeof(logtrexp), M::AbstractArray{T}) where T<:Real
    e, v = symeigen(M)
    y = logsumexp(e)
    function logtrexp_pullback(ȳ)
        Λ = map(x->exp(x-y), e) 
        ∂y_∂M = v * diagm(Λ) * v' |> transpose
        M̄ = ȳ * ∂y_∂M
        return ChainRules.NoTangent(), M̄
    end
    return y, logtrexp_pullback
end

#end #module rrule
