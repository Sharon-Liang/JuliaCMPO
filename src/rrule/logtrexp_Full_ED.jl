#module rrule
"""
    rrule for logtrexp(tM) function where `typeof(M) <: AbstractArray`
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


"""
    rrule for logtrexp(tM) function where `typeof(M)=CMPSMatrix`, 
    `∂y_∂M` also return a `CMPSMatrix`. `exp(tM)` is constructed 
    explicitly if no `estimator` is assigned.
"""
function ChainRules.rrule(::typeof(logtrexp), t::Real, M::CMPSMatrix{Ts,T,S,U}) where {Ts,T,S,U}
    @unpack ψl, ψr = M
    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)

    M = Matrix(M) |> symmetrize
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

        if U <: AbstractMatrix
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
        ψ̄l = CMPS_generate(Q̄l, R̄l)
        ψ̄r = CMPS_generate(Q̄r, R̄r)
        M̄ = CMPSMatrix(ψ̄l, ψ̄r)
        return ChainRules.NoTangent(), t̄, M̄
    end
    return y, logtrexp_pullback
end


"""
    rrule for logtrexp(tM) function where `typeof(M)=CMPSMatrix`, 
    `∂y_∂M` also return a `CMPSMatrix`. When `trace_estimator = Full_ED`, 
    `exp(tM)` is not explicitly constructed. 
    (Note there is no memory saving in this case because one has to 
    save all eigen vectors of the CMPSMatrix.)
"""
function ChainRules.rrule(::typeof(logtrexp), 
                          t::Real, M::CMPSMatrix{Ts,T,S,U}, 
                          trace_estimator::TraceEstimator{Tf, To}
                          ) where {Ts,T,S,U,Tf<:typeof(Full_ED),To}
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

        χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)
        if processor == CPU
            Onel = ones(T, χl, χl)
            Oner = ones(T, χr, χr)
        else
            Onel = CUDA.ones(T, χl, χl)
            Oner = CUDA.ones(T, χr, χr)
        end

        ∂y_∂Ql = -t * ein"n,kbn,kan,kk -> ab"(Λ, conj(vecs), vecs, Oner)
        ∂y_∂Qr = -t * ein"n,bkn,akn,kk -> ab"(Λ, conj(vecs), vecs, Onel)
        if U <: AbstractMatrix
            ∂y_∂Rl = -t * ein"n,kl,lbn,kan,kl -> ab"(Λ, ψr.R, conj(vecs), vecs, Oner)
            ∂y_∂Rr = -t * ein"n,kl,bln,akn,kl -> ab"(Λ, ψl.R, conj(vecs), vecs, Onel)
        else
            if processor == CPU
                Onel = ones(T, χl, χl, size(ψl.R,3))
                Oner = ones(T, χr, χr, size(ψr.R,3))
            else
                Onel = CUDA.ones(T, χl, χl, size(ψl.R,3))
                Oner = CUDA.ones(T, χr, χr, size(ψr.R,3))
            end
            ∂y_∂Rl = -t * ein"n,klm,lbn,kan,klm -> abm"(Λ, ψr.R, conj(vecs), vecs, Oner)
            ∂y_∂Rr = -t * ein"n,klm,bln,akn,klm -> abm"(Λ, ψl.R, conj(vecs), vecs, Onel)
        end
        Q̄l = ȳ * ∂y_∂Ql
        R̄l = ȳ * ∂y_∂Rl
        Q̄r = ȳ * ∂y_∂Qr
        R̄r = ȳ * ∂y_∂Rr
        ψ̄l = CMPS_generate(Q̄l, R̄l)
        ψ̄r = CMPS_generate(Q̄r, R̄r)
        M̄ = CMPSMatrix(ψ̄l, ψ̄r)
        return ChainRules.NoTangent(), t̄, M̄, ChainRules.NoTangent()
    end
    return y, logtrexp_pullback
end


#end #module rrule
