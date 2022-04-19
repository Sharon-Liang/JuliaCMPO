#module cMPSOperations
"""
    Logarithm of the overlap of two CMPS:
        log_overlap = ln(⟨ψl|ψr⟩)
"""
function log_overlap(sl::T, sr::T, β::Real) where T<:AbstractCMPS
    K = sl * sr
    return logtrexp(-β * K)
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
function normalize(s::AbstractCMPS, β::Real)
    λ = log_overlap(s, s, β)/β
    eye = λ/2 * Matrix{Float64}(I,size(s.Q))
    Q = s.Q - convert(typeof(s.Q), eye)
    return CMPS_generate(Q, s.R)
end


"""
    transpose of a CMPO
"""
function transpose(o::AbstractCMPO)
    Q = o.Q
    R = o.L
    L = o.R
    length(size(o.P)) == 2 ? P = ein"ij->ij"(o.P) : P = ein"ijkl->ijlk"(o.P)
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
function project(s::T, u::AbstractMatrix) where T<:AbstractCMPS
    u = convert(typeof(s.Q), u)
    au = convert(typeof(s.Q), u')
    Q = au * s.Q * u
    if length(size(s.R)) == 2
        R = au * s.R * u
    else
        R = ein"(ip,pql),qj -> ijl"(au, s.R, u)
    end
    return CMPS_generate(Q, R)
end


function project(o::AbstractCMPO, u::AbstractMatrix) where T<:AbstractCMPO
    u = convert(typeof(o.Q), u)
    au = convert(typeof(o.Q), u')
    Q = au * o.Q * u
    if length(size(o.R)) == 2
        R = au * o.R * u
        L = au * o.R * u
        P = au * o.P * u
    else
        R = ein"(ip,pql),qj -> ijl"(au, o.R, u)
        L = ein"(ip,pql),qj -> ijl"(au, o.L, u)
        P = ein"(ip,pqlm),qj -> ijlm"(au, o.P, u)
    end
    return CMPO_generate(Q,R,L,P)
end


""" multiply a matrix to the left/right of the cMPS/cMPO
    --        --   --      --
    | 1 0 ... 0|   | I + ϵQ |
    | 0        |   |        |
    | :        |   |   |    |
    | :    W   |   |  √ϵR   |
    | 0        |   |   |    |
    --        --   --      --
"""
function *(Ut::AbstractMatrix, s::AbstractCMPS)
    Ut = convert(typeof(s.Q), Ut)
    if length(size(s.R)) == 2 
        R = CUDA.@allowscalar Ut[1,1] * s.R
    else
        R = ein"mn, ijn -> ijm"(Ut, s.R)
    end
    return CMPS_generate(s.Q, R)
end

function *(Ut::AbstractMatrix, o::AbstractCMPO)
    Ut = convert(typeof(o.Q), Ut)
    if length(size(o.R)) == 2
        R = CUDA.@allowscalar Ut[1,1] * o.R
        P = CUDA.@allowscalar Ut[1,1] * o.P
    else
        R = ein"mk, ijk -> ijm"(Ut, o.R)
        P = ein"mk, ijkn -> ijmn"(Ut, o.P)
    end
    return CMPO_generate(o.Q, R, o.L, P)
end

function *(o::AbstractCMPO, Ut::AbstractMatrix)
    Ut = convert(typeof(o.Q), Ut)
    if length(size(o.L)) == 2
        L = CUDA.@allowscalar Ut[1,1] * o.L
        P = CUDA.@allowscalar Ut[1,1] * o.P
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
    _, u = symeigen(s.Q)
    return project(s, u)
end


#end    #module cMPSOperations
