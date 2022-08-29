#module rrule
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


"""
    ChainRules.rrule(::typeof(logtrexp), t::Real, M::AbstractArray)

Rrule for ``ln[Tr(e^{tM})]`` function where `typeof(M) <: CMPSMatrix`.
`∂y_∂M` is also of type `CMPSMatrix`. `exp(tM)` is constructed explicitly 
if no `trace_estimator` is assigned.
"""
function ChainRules.rrule(::typeof(logtrexp), t::Real, M::CMPSMatrix)
    @unpack ψl, ψr = M
    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)

    M = tomatrix(M) |> symmetrize
    vals, vecs = eigensolver(M)
    y = logsumexp(t*vals)

    function logtrexp_pullback(ȳ)
        ∂y_∂t = map(e -> e * exp(t*e - y), vals) |> sum
        t̄ = ȳ * ∂y_∂t

        Λ = map(e -> exp(t*e - y), vals) 
        expm = vecs * diagm(Λ) * vecs'
        expm = reshape(expm, χr, χl, χr, χl)
        ∂y_∂Ql = -t * ein"kbka -> ab"(expm)
        ∂y_∂Qr = -t * ein"bkak -> ab"(expm)

        if typeof(ψl.R) <: AbstractMatrix
            ∂y_∂Rl = -t * ein"kl, lbka -> ab"(ψr.R, expm)
            ∂y_∂Rr = -t * ein"kl, blak -> ab"(ψl.R, expm)
        else
            ∂y_∂Rl = -t * ein"klm, lbka -> abm"(ψr.R, expm)
            ∂y_∂Rr = -t * ein"klm, blak -> abm"(ψl.R, expm)
        end
        Q̄l = ȳ * ∂y_∂Ql
        R̄l = ȳ * ∂y_∂Rl
        Q̄r = ȳ * ∂y_∂Qr
        R̄r = ȳ * ∂y_∂Rr
        ψ̄l = cmps_generate(Q̄l, R̄l)
        ψ̄r = cmps_generate(Q̄r, R̄r)
        M̄ = CMPSMatrix(ψ̄l, ψ̄r)
        return ChainRules.NoTangent(), t̄, M̄
    end
    return y, logtrexp_pullback
end


"""
ChainRules.rrule(::typeof(logtrexp), t::Real, M::AbstractArray, trace_estimator::TraceEstimator{<:typeof(full_ed)})

Rrule for ``ln[Tr(e^{tM})]`` function where `typeof(M) <: CMPSMatrix`.
`∂y_∂M` is also of type `CMPSMatrix`. When `trace_estimator = full_ed_ftrace`, 
`exp(tM)` is not explicitly constructed. Note that there is no memory saving 
in this case because one has to save all eigen vectors of the `CMPSMatrix.
"""
function ChainRules.rrule(::typeof(logtrexp), t::Real, M::CMPSMatrix, 
                          trace_estimator::TraceEstimator{<:typeof(full_ed_ftrace), To}) where {To}
    @unpack ψl, ψr = M
    @unpack processor = trace_estimator.options
    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)
    vals, vecs = eigensolver(M)
    y = logsumexp(t*vals)

    Ns = size(M, 1)
    function logtrexp_pullback(ȳ)
        ∂y_∂t = map(e -> e * exp(t*e - y), vals) |> sum
        t̄ = ȳ * ∂y_∂t
        
        Λ = map(e -> exp(t*e - y), vals) 
        vecs = reshape(vecs, χr, χl, Ns)

        ∂y_∂Ql = -t * ein"n,kbn,kan -> ab"(Λ, conj(vecs), vecs)
        ∂y_∂Qr = -t * ein"n,bkn,akn -> ab"(Λ, conj(vecs), vecs)
        if typeof(ψl.R) <: AbstractMatrix
            ∂y_∂Rl = -t * ein"n,kl,lbn,kan -> ab"(Λ, ψr.R, conj(vecs), vecs)
            ∂y_∂Rr = -t * ein"n,kl,bln,akn -> ab"(Λ, ψl.R, conj(vecs), vecs)
        else
            ∂y_∂Rl = -t * ein"n,klm,lbn,kan -> abm"(Λ, ψr.R, conj(vecs), vecs)
            ∂y_∂Rr = -t * ein"n,klm,bln,akn -> abm"(Λ, ψl.R, conj(vecs), vecs)
        end
        Q̄l = ȳ * ∂y_∂Ql
        R̄l = ȳ * ∂y_∂Rl
        Q̄r = ȳ * ∂y_∂Qr
        R̄r = ȳ * ∂y_∂Rr
        ψ̄l = cmps_generate(Q̄l, R̄l)
        ψ̄r = cmps_generate(Q̄r, R̄r)
        M̄ = CMPSMatrix(ψ̄l, ψ̄r)
        return ChainRules.NoTangent(), t̄, M̄, ChainRules.NoTangent()
    end
    return y, logtrexp_pullback
end


#end #module rrule
