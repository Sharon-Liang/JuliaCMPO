#module cMPSOperations
import Base: ==, ≈, *
"""
    tomatrix(A::CMPSMatrix)

Generate the matrix form of `A`.
"""
function tomatrix(A::CMPSMatrix) 
    @unpack ψl, ψr = A
    li, ri = oneunit(ψl.Q), oneunit(ψr.Q)
    K = li ⊗ ψr.Q + ψl.Q ⊗ ri + ψl.R ⊗ ψr.R
    return -K
end


"""
    log_overlap(ψl, ψr, β[, trace_estimator])

Calculate ``ln(⟨ψl|ψr⟩)``
"""
function log_overlap(ψl::AbstractCMPS, ψr::AbstractCMPS, β::Real, trace_estimator = nothing)
    K = CMPSMatrix(ψl, ψr)
    return logtrexp(-β, K, trace_estimator)
end


"""
    norm(s::AbstractCMPS, β[, trace_estimator])

Calculate the norm of a CMPS `s`.
"""
function LinearAlgebra.norm(s::AbstractCMPS, β::Real, trace_estimator = nothing)
    λ = log_overlap(s, s, β, trace_estimator) 
    return exp(λ/2)
end


"""
    normalize(s::AbstractCMPS, β[, trace_estimator])

Normalize a CMPS `s`, i.e. make ``⟨ψ|ψ⟩ = 1``.
"""
function LinearAlgebra.normalize(s::AbstractCMPS, β::Real, trace_estimator = nothing)
    λ = log_overlap(s, s, β, trace_estimator)/β
    Q = s.Q - λ/2 * oneunit(s.Q)
    return CMPS_generate(Q, s.R)
end


"""
    logfidelity(ψ, ψ0, β[, trace_estimator])

Calculate the logarithm of fidelity function ``\\mathcal{F} = \\frac{⟨ψ|ψ0⟩}{√⟨ψ|ψ⟩}``
between two CMPSs `ψ` and `ψ0`.
"""
function logfidelity(ψ::AbstractCMPS, ψ0::AbstractCMPS, β::Real, trace_estimator = nothing)
    return log_overlap(ψ, ψ0, β, trace_estimator) - 0.5*log_overlap(ψ, ψ, β, trace_estimator)
end


"""
    fidelity(ψ, ψ0, β[, trace_estimator; ifnormalize::Bool=false])

Calculate the fidelity function ``\\mathcal{F} = \\frac{⟨ψ|ψ0⟩}{√⟨ψ|ψ⟩}``
between two CMPSs `ψ` and `ψ0`.
"""
function fidelity(ψ::AbstractCMPS, ψ0::AbstractCMPS, β::Real, trace_estimator = nothing; 
                  ifnormalize::Bool = false)
    if ifnormalize
        ψ = normalize(ψ, β, trace_estimator)
        ψ0 = normalize(ψ0, β, trace_estimator)
    end
    return logfidelity(ψ, ψ0, β, trace_estimator) |> exp
end


"""
    transpose(o::AbstractCMPO)
    
Calculate the transpose of a CMPO.
"""
function Base.transpose(o::AbstractCMPO)
    Q = o.Q
    R = o.L
    L = o.R
    typeof(o.P) <: AbstractMatrix ? P = ein"ij->ij"(o.P) : P = ein"ijkl->ijlk"(o.P)
    return CMPO_generate(Q,R,L,P)
end


"""
    adjoint(o::AbstractCMPO)
    
Calculate the hermitian conjuate (adjoint) of a CMPO.
"""
function Base.adjoint(o::AbstractCMPO)
    ot = transpose(o)
    Q = conj.(ot.Q)
    R = conj.(ot.L)
    L = conj.(ot.R)
    P = conj.(ot.P)
    return CMPO_generate(Q,R,L,P)
end


"""
    ==(a::AbstractCTensor, b::AbstractCTensor)

Determine if two AbstractCTensors are equal.
"""
function ==(a::AbstractCMPO, b::AbstractCMPO) 
    return ==(a.Q, b.Q) && 
            ==(a.R, b.R) &&
            ==(a.L, b.L) &&
            ==(a.P, b.P) 
end

function ==(a::AbstractCMPS, b::AbstractCMPS)
    return ==(a.Q, b.Q) && ==(a.R, b.R)
end


"""
    ==(a::CMPSMatrix, b::CMPSMatrix)

    Determine if two CMPSMatrices are equal.
"""
function ==(a::CMPSMatrix, b::CMPSMatrix) 
    return ==(a.ψl, b.ψl) && ==(a.ψr, b.ψr)
end


"""
    ≈(a::AbstractCTensor, b::AbstractCTensor)
"""
function ≈(a::AbstractCMPO, b::AbstractCMPO; kwargs...) 
    return ≈(a.Q, b.Q; kwargs...) && 
            ≈(a.R, b.R; kwargs...) &&
            ≈(a.L, b.L; kwargs...) &&
            ≈(a.P, b.P; kwargs...) 
end

function ≈(a::AbstractCMPS, b::AbstractCMPS; kwargs...) 
    return ≈(a.Q, b.Q; kwargs...) && 
           ≈(a.R, b.R; kwargs...)
end


"""
    ≈(a::CMPSMatrix, b::CMPSMatrix) 
"""
function ≈(a::CMPSMatrix, b::CMPSMatrix; kwargs...)
    return ≈(a.ψl, b.ψl; kwargs...) && ≈(a.ψr, b.ψr; kwargs...)
end


"""
    ishermitian(o::AbstractCMPO)

Determine if a CMPO is hermitian.
"""
LinearAlgebra.ishermitian(o::AbstractCMPO) = isequal(o, adjoint(o))


"""
    project(s::AbstractCTensor, u::AbstractMatrix)

Perform a unitary transformation in the imaginary-time direction of a CMPS
(i.e. On the physical bond of CMPO and virtual bond of CMPS).
This is a guage transformation if `u` is a square matrix.
"""
function project(s::AbstractCMPS, u::AbstractMatrix) 
    Q = u' * s.Q * u
    if typeof(s.R) <: AbstractMatrix
        R = u' * s.R * u
    else
        R = ein"(ip,pql),qj -> ijl"(u', s.R, u)
    end
    return CMPS_generate(Q, R)
end

function project(o::AbstractCMPO, u::AbstractMatrix)
    Q = u' * o.Q * u
    if typeof(o.R) <: AbstractMatrix
        R = u' * o.R * u
        L = u' * o.R * u
        P = u' * o.P * u
    else
        R = ein"(ip,pql),qj -> ijl"(u', o.R, u)
        L = ein"(ip,pql),qj -> ijl"(u', o.L, u)
        P = ein"(ip,pqlm),qj -> ijlm"(u', o.P, u)
    end
    return CMPO_generate(Q,R,L,P)
end


""" 
    *(Ut::AbstractMatrix, s::AbstractCMPS)
    *(s::AbstractCMPS, Ut::AbstractMatrix)

multiply a matrix to the left/right of a CMPS
    --        --   --      --
    | 1 0 ... 0|   | I + ϵQ |
    | 0        |   |        |
    | :        |   |   |    |
    | :    W   |   |  √ϵR   |
    | 0        |   |   |    |
    --        --   --      --
"""
function *(Ut::AbstractMatrix, s::AbstractCMPS)
    if typeof(s.R) <: AbstractMatrix
        R = sum(Ut) * s.R
    else
        R = ein"mn, ijn -> ijm"(Ut, s.R)
    end
    return CMPS_generate(s.Q, R)
end


function *(Ut::AbstractMatrix, o::AbstractCMPO)
    if typeof(o.R) <: AbstractMatrix
        R = sum(Ut) * o.R
        P = sum(Ut) * o.P
    else
        R = ein"mk, ijk -> ijm"(Ut, o.R)
        P = ein"mk, ijkn -> ijmn"(Ut, o.P)
    end
    return CMPO_generate(o.Q, R, o.L, P)
end


function *(o::AbstractCMPO, Ut::AbstractMatrix)
    if typeof(o.R) <: AbstractMatrix
        L = sum(Ut) * o.L
        P = sum(Ut) * o.P
    else
        L = ein"ijk, km -> ijm"(o.L, Ut)
        P = ein"ijnk, km -> ijnm"(o.P, Ut) 
    end
    return CMPO_generate(o.Q, o.R, L, P)
end


""" 
    diagQ(s::AbstractCMPS)

Transform the CMPS `s` to the gauge where `s.Q` is a diagonalized matrix.
"""
function diagQ(s::AbstractCMPS)
    Q = symmetrize(s.Q)
    _, u = eigensolver(Q)
    return project(s, u)
end


#end    #module cMPSOperations
