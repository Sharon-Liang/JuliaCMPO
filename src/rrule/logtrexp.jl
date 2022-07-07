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


function ChainRules.rrule(::typeof(logtrexp), t::Real, M::CMPSMatrix{T, S, U}) where {T,S,U}
    @unpack ψl, ψr = M
    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)

    M = Matrix(M)
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


#=
"""
    rrule for logtrexp(tM) function where `typeof(M)=CMPSMatrix`, 
    `∂y_∂M` also return a `CMPSMatrix`.
"""
function ChainRules.rrule(::typeof(logtrexp), t::Real, M::T, method::typeof(Full_ED)) where {T<:CMPSMatrix}
    sign(t) == 1 ? which = :LR : which=:SR
    e0, _, _ = eigsolve(M, size(M,1), 1, which, ishermitian = true)
    e0 = e0[1]
    expr = e -> exp(t * (e - e0))
    res = method(M, expr)[1]
    return log(res) + t*e0
        e0, _, _ = eigsolve(M, size(M,1), 1, :SR, ishermitian = true)
        e0 = e0[1]
        exprZ = e -> exp(-β*(e - e0))
        return method(M, exprZ)[1]
    end
    
    
    
    vals, vecs = eigensolver(M)
    y = logsumexp(t*vals)
    function logtrexp_pullback(ȳ)
        ∂y_∂t = map(e -> e * exp(t*e - y), vals) |> sum
        t̄ = ȳ * ∂y_∂t

        Λ = map(e -> exp(t*e - y), vals) 
        ∂y_∂M = vecs * diagm(Λ) * vecs' |> transpose
        M̄ = ȳ * t * ∂y_∂M
        return ChainRules.NoTangent(), t̄, M̄, ChainRules.NoTangent()
    end
    return y, logtrexp_pullback
end
=#
#end #module rrule
