#module cMPSOperations
"""
    `Matrix(A::CMPSMatrix{T, S, U})`: return the matrix form of A
"""
function Matrix(A::CMPSMatrix{Ts,T,S,U}) where {Ts,T,S,U}
    @unpack ψl, ψr = A
    li = convert(S, Matrix{T}(I,size(ψl.Q)))
    ri = convert(S, Matrix{T}(I,size(ψr.Q)))
    K = li ⊗ ψr.Q + ψl.Q ⊗ ri + ψl.R ⊗ ψr.R
    return -K
end


"""
    Logarithm of the overlap of two CMPS:
        log_overlap = ln(⟨ψl|ψr⟩)
"""
function log_overlap(ψl::T, ψr::T, β::Real) where T<:AbstractCMPS
    K = ψl * ψr
    return logtrexp(-β, K)
end


"""
    Norm of a CMPS
"""
function norm(s::AbstractCMPS, β::Real)
    λ = log_overlap(s, s, β) 
    return exp(λ/2)
end


"""
    Normalize a CMPS, i.e. ⟨ψ|ψ⟩ = 1
"""
function normalize(s::AbstractCMPS{T, S, U}, β::Real) where {T,S,U}
    λ = log_overlap(s, s, β)/β
    eye = λ/2 * Matrix{Float64}(I,size(s.Q))
    Q = s.Q - convert(S, eye)
    return CMPS_generate(Q, s.R)
end


"""
    Fidelity between the target cMPS |ψ⟩ and the origional cMPS T|r⟩:
        fidelity = ⟨ψ|T|r⟩/√(⟨ψ|ψ⟩)
        logfidelity(ψ, ψ0) = ln(Fd)
"""
logfidelity(ψ::T, ψ0::T, β::Real) where T<:AbstractCMPS = log_overlap(ψ, ψ0, β) - 0.5*log_overlap(ψ, ψ, β)

function fidelity(ψ::T, ψ0::T, β::Real; Normalize::Bool = false) where T<:AbstractCMPS
    if Normalize
        ψ = normalize(ψ, β)
        ψ0 = normalize(ψ0, β)
    end
    return logfidelity(ψ, ψ0, β) |> exp
end


"""
    transpose of a CMPO
"""
function transpose(o::AbstractCMPO{T, S, U, V}) where {T,S,U,V}
    Q = o.Q
    R = o.L
    L = o.R
    U <: AbstractMatrix ? P = ein"ij->ij"(o.P) : P = ein"ijkl->ijlk"(o.P)
    return CMPO_generate(Q,R,L,P)
end


"""
    Hermitian conjuate (adjoint) of a CMPO
"""
function adjoint(o::AbstractCMPO)
    ot = transpose(o)
    Q = conj.(ot.Q)
    R = conj.(ot.L)
    L = conj.(ot.R)
    P = conj.(ot.P)
    return CMPO_generate(Q,R,L,P)
end


"""
    ==(a::T, b::T) where T<:AbstractCTensor
"""
function ==(a::T, b::T) where T<:AbstractCMPO
    return ==(a.Q, b.Q) && 
            ==(a.R, b.R) &&
            ==(a.L, b.L) &&
            ==(a.P, b.P) 
end

function ==(a::T, b::T) where T<:AbstractCMPS
    return ==(a.Q, b.Q) && ==(a.R, b.R)
end


"""
    ≈(a::T, b::T) where T<:AbstractCTensor
"""
function ≈(a::T, b::T; kwargs...) where T<:AbstractCMPO
    return ≈(a.Q, b.Q; kwargs...) && 
            ≈(a.R, b.R; kwargs...) &&
            ≈(a.L, b.L; kwargs...) &&
            ≈(a.P, b.P; kwargs...) 
end

function ≈(a::T, b::T; kwargs...) where T<:AbstractCMPS
    return ≈(a.Q, b.Q; kwargs...) && ≈(a.R, b.R)
end


"""
    Determine if a CMPO is hermitian
"""
ishermitian(o::AbstractCMPO) = isequal(o, adjoint(o))


"""
    Perform a unitary transformation in the imaginary-time direction
    (i.e. On the physical bond of cMPO and virtual bond of cMPS)
    if `U` is a square matrix, this is a guage transformation
"""
function project(s::AbstractCMPS{T,S,U}, u::AbstractMatrix) where {T,S,U}
    Q = u' * s.Q * u
    if U <: AbstractMatrix
        R = u' * s.R * u
    else
        R = ein"(ip,pql),qj -> ijl"(u', s.R, u)
    end
    return CMPS_generate(Q, R)
end

function project(o::AbstractCMPO{T, S, U, V}, u::AbstractMatrix) where {T,S,U,V}
    Q = u' * o.Q * u
    if U <: AbstractMatrix
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


""" multiply a matrix to the left/right of the cMPS
    --        --   --      --
    | 1 0 ... 0|   | I + ϵQ |
    | 0        |   |        |
    | :        |   |   |    |
    | :    W   |   |  √ϵR   |
    | 0        |   |   |    |
    --        --   --      --
"""
function *(Ut::AbstractMatrix, s::AbstractCMPS{T,S,U}) where {T,S,U}
    if U <: AbstractMatrix
        R = sum(Ut) * s.R
    else
        R = ein"mn, ijn -> ijm"(Ut, s.R)
    end
    return CMPS_generate(s.Q, R)
end

function *(Ut::AbstractMatrix, o::AbstractCMPO{T, S, U, V}) where {T,S,U,V}
    if U <: AbstractMatrix
        R = sum(Ut) * o.R
        P = sum(Ut) * o.P
    else
        R = ein"mk, ijk -> ijm"(Ut, o.R)
        P = ein"mk, ijkn -> ijmn"(Ut, o.P)
    end
    return CMPO_generate(o.Q, R, o.L, P)
end

function *(o::AbstractCMPO{T, S, U, V}, Ut::AbstractMatrix) where {T,S,U,V}
    if U <: AbstractMatrix
        L = sum(Ut) * o.L
        P = sum(Ut) * o.P
    else
        L = ein"ijk, km -> ijm"(o.L, Ut)
        P = ein"ijnk, km -> ijnm"(o.P, Ut) 
    end
    return CMPO_generate(o.Q, o.R, L, P)
end

""" 
    transform the CMPS to the gauge where Q is a diagonalized matrix
"""
function diagQ(s::AbstractCMPS)
    Q = symmetrize(s.Q)
    _, u = eigensolver(Q)
    return project(s, u)
end


#end    #module cMPSOperations
