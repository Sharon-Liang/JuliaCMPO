#=
### *Multiplications* : *Otimes* 
=#

"""
    ⊗(A, B)

Multiplications between arrays in CMPS and CMPO data structures
"""
function ⊗(A::AbstractMatrix, B::AbstractMatrix)
    r1, c1 = size(A)
    r2, c2 = size(B)
    return reshape(ein"ij,kl->kilj"(A, B), r1*r2, c1*c2)
end


function ⊗(A::AbstractMatrix, B::AbstractArray{T,3} where T) 
    r1, c1 = size(A)
    r2, c2, d = size(B)
    return reshape(ein"ij,klm->kiljm"(A, B), r1*r2, c1*c2, d)
end


function ⊗(A::AbstractArray{T,3} where T, B::AbstractMatrix)
    r1, c1, d = size(A)
    r2, c2 = size(B)
    return reshape(ein"ijm,kl->kiljm"(A, B), r1*r2, c1*c2, d)
end


function ⊗(A::AbstractArray{T,3} where T, B::AbstractArray{T,3} where T)
    r1, c1, _ = size(A)
    r2, c2, _ = size(B)
    return reshape(ein"ijm,klm->kilj"(A, B), r1*r2, c1*c2)
end


function ⊗(A::AbstractArray{T,4} where T, B::AbstractArray{T,3} where T)
    r1, c1, d, _ = size(A)
    r2, c2, _ = size(B)
    return reshape(ein"ijnm,klm->kiljn"(A, B), r1*r2, c1*c2, d)    
end


function ⊗(A::AbstractArray{T,3} where T, B::AbstractArray{T,4} where T)
    r1, c1, _ = size(A)
    r2, c2, _, f = size(B)
    return reshape(ein"ijm,klmn->kiljn"(A, B), r1*r2, c1*c2, f)    
end


function ⊗(A::AbstractArray{T,4} where T, B::AbstractArray{T,4} where T)
    r1, c1, d, _ = size(A)
    r2, c2, _, f = size(B)
    return reshape(ein"ijpm,klmq->kiljpq"(A, B), r1*r2, c1*c2, d, f)    
end


#=
### *Multiplications* : *Products between CMPS and CMPO* 
=#
import Base:*

"""
    *(a::Union{CMPS, CMPO}, b::Union{CMPS, CMPO})

Products between `CMPS` and `CMPO`.
"""
function *(sl::CMPS, sr::CMPS)
    li, ri = oneunit(sl.Q), oneunit(sr.Q)
    K = li ⊗ sr.Q + sl.Q ⊗ ri + sl.R ⊗ sr.R
    return -K
end


function *(o::CMPO, s::CMPS)
    oi, si = oneunit(o.Q), oneunit(s.Q)
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return CMPS(Q, R)
end


function *(s::CMPS, o::CMPO)
    oi, si = oneunit(o.Q), oneunit(s.Q)
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return CMPS(Q, R)
end


function *(ol::CMPO, or::CMPO)
    li, ri = oneunit(ol.Q), oneunit(or.Q)
    Q = ol.Q ⊗ ri + li ⊗ or.Q + ol.L ⊗ or.R
    L = li ⊗ or.L + ol.L ⊗ or.P
    R = ol.R ⊗ ri + ol.P ⊗ or.R
    P = ol.P ⊗ or.P
    return CMPO(Q,R,L,P)
end


#=
### *Properties* : *CMPS* 
=#
"""
    log_overlap(ψl, ψr, β)

Calculate ``ln(⟨ψl|ψr⟩)``.
"""
function log_overlap(ψl::CMPS, ψr::CMPS, β::Real)
    K = ψl * ψr
    return logtrexp(-β, K)
end


"""
    norm(s::CMPS, β)

Calculate the norm of a cMPS.
"""
function LinearAlgebra.norm(s::CMPS, β::Real)
    λ = log_overlap(s, s, β) 
    return exp(λ/2)
end


"""
    normalize(s::CMPS, β)

Normalize a cMPS, i.e. to make ``⟨ψ|ψ⟩ = 1``.
"""
function LinearAlgebra.normalize(s::CMPS, β::Real)
    λ = log_overlap(s, s, β)/β
    Q = s.Q - λ/2 * oneunit(s.Q)
    return CMPS(Q, s.R)
end


"""
    logfidelity(ψ, ψ₀, β[, tonormalize::Bool = true])

Calculate the logarithm of fidelity function ``\\mathcal{F} = \\frac{⟨ψ|ψ₀⟩}{√⟨ψ|ψ⟩}`` between two cMPS `ψ` and `ψ₀`.
"""
function logfidelity(ψ::CMPS, ψ₀::CMPS, β::Real, tonormalize::Bool = true) 
    if tonormalize
        ψ = normalize(ψ, β)
        ψ₀ = normalize(ψ₀, β)
    end
    return log_overlap(ψ, ψ₀, β) - 0.5*log_overlap(ψ, ψ, β)
end


"""
    fidelity(ψ, ψ₀, β[; tonormalize::Bool = false])

Calculate the fidelity function ``\\mathcal{F} = \\frac{⟨ψ|ψ₀⟩}{√⟨ψ|ψ⟩}`` between two cMPS `ψ` and `ψ₀`.
"""
fidelity(ψ::CMPS, ψ₀::CMPS, β::Real, tonormalize::Bool = true) = logfidelity(ψ, ψ₀, β, tonormalize) |> exp


"""
    project(s::AbstractCTensor, u::AbstractMatrix)

Perform a unitary transformation in the imaginary-time direction of a CMPS. This is a guage transformation if `u` is a square matrix.
"""
function project(s::CMPS, u::AbstractMatrix) 
    Q = u' * s.Q * u
    if typeof(s.R) <: AbstractMatrix
        R = u' * s.R * u
    else
        R = ein"(ip,pql),qj -> ijl"(u', s.R, u)
    end
    return CMPS(Q, R)
end


""" 
    diagQ(s::CMPS)

Transform the CMPS `s` to the gauge where `s.Q` is a diagonalized matrix.
"""
function diagQ(s::CMPS)
    Q = symmetrize(s.Q)
    _, u = eigensolver(Q)
    return project(s, u)
end


#=
### *Properties* : *CMPO* 
=#
"""
    transpose(o::CMPO)
    
Calculate the transpose of a cMPO.
"""
function Base.transpose(o::CMPO)
    Q = o.Q
    R = o.L
    L = o.R
    typeof(o.P) <: AbstractMatrix ? P = ein"ij->ij"(o.P) : P = ein"ijkl->ijlk"(o.P)
    return CMPO(Q,R,L,P)
end


"""
    adjoint(o::CMPO)
    
Calculate the hermitian conjuate (adjoint) of a CMPO.
"""
function Base.adjoint(o::CMPO)
    ot = transpose(o)
    Q = conj.(ot.Q)
    R = conj.(ot.L)
    L = conj.(ot.R)
    P = conj.(ot.P)
    return CMPO(Q,R,L,P)
end


#=
### *Compare*
=#
import Base: ==, ≈

"""
    ==(a, b)

Determine if two cMPS/cMPO are equal.
"""
function ==(a::CMPO, b::CMPO) 
    return ==(a.Q, b.Q) && 
            ==(a.R, b.R) &&
            ==(a.L, b.L) &&
            ==(a.P, b.P) 
end

function ==(a::CMPS, b::CMPS)
    return ==(a.Q, b.Q) && ==(a.R, b.R)
end


"""
    ≈(a, b[; kwargs...])

Determine if two cMPS/cMPO are approximate.
"""
function ≈(a::CMPO, b::CMPO; kwargs...) 
    return ≈(a.Q, b.Q; kwargs...) && 
            ≈(a.R, b.R; kwargs...) &&
            ≈(a.L, b.L; kwargs...) &&
            ≈(a.P, b.P; kwargs...) 
end

function ≈(a::CMPS, b::CMPS; kwargs...) 
    return ≈(a.Q, b.Q; kwargs...) && 
           ≈(a.R, b.R; kwargs...)
end


"""
    ishermitian(o::CMPO)

Determine if a cMPO is hermitian.
"""
LinearAlgebra.ishermitian(o::CMPO) = isequal(o, adjoint(o))


