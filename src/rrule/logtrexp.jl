#module rrule
Zygote.@adjoint Array(x::CuArray) = Array(x), dy->(CuArray(dy),)

"""
    rrule for logtrexp(tM) function
"""    
function ChainRules.rrule(::typeof(logtrexp), t::Real, M::T) where {T<:AbstractArray}
    M = M |> symmetrize
    vals, vecs = eigensolver(M)
    y = logsumexp(t*vals)
    function logtrexp_pullback(ȳ)
        ∂y_∂t = map(e -> e * exp(t*e - y), vals) |> sum
        t̄ = ȳ * ∂y_∂t

        Λ = map(e -> exp(t*e - y), vals) 
        ∂y_∂M = vecs * diagm(Λ) * vecs' |> transpose
        M̄ = ȳ * t * ∂y_∂M
        return ChainRules.NoTangent(), t̄, M̄
    end
    return y, logtrexp_pullback
end


#end #module rrule
