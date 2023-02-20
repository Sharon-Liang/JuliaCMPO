ChainRulesCore.@non_differentiable oneunit(::Any...)

"""
    ChainRules.rrule(::typeof(logtrexp), t::Real, M::AbstractArray)

Rrule for ``ln[Tr(e^{tM})]`` function where `typeof(M) <: AbstractArray`.
"""    
function ChainRules.rrule(::typeof(logtrexp), t::Real, M::AbstractArray)
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