#module cMPSOperations
"""
    Logarithm of the overlap of two CMPS:
        log_overlap = ln(⟨ψl|ψr⟩)
"""
function log_overlap(sl::CMPS, sr::CMPS, β::Real)
    K = sl * sr
    return logtrexp(-β * K)
end


"""
    Norm of a CMPS
"""
function norm(s::CMPS, β::Real)
    λ = log_overlap(s, s, β)
    return exp(β*λ/2)
end


"""
    Normalize a CMPS, i.e. ⟨ψ|ψ⟩ = 1
"""
function normalize(s::CMPS, β::Real)
    λ = log_overlap(s, s, β)/β
    eye = λ/2 * Matrix{Float64}(I,size(s.Q))
    Q = s.Q - eye
    return CMPS(Q, s.R)
end


"""
    transpose of a CMPO
"""
function transpose(o::CMPO)
    Q = o.Q
    R = o.L
    L = o.R
    length(size(o.P)) == 2 ? P = ein"ij->ij"(o.P) : P = ein"ijkl->ijlk"(o.P)
    return CMPO(Q,R,L,P)
end


"""
    Hermitian conjuate (adjoint) of a CMPO
"""
function adjoint(o::CMPO)
    ot = transpose(o)
    Q = conj.(ot.Q)
    R = conj.(ot.L)
    L = conj.(ot.R)
    P = conj.(ot.P)
    return CMPO(Q,R,L,P)
end


"""
    isequal(a::CMPO, b::CMPO)
"""
function isequal(a::CMPO, b::CMPO)
    return isequal(a.Q, b.Q) && 
            isequal(a.R, b.R) &&
            isequal(a.L, b.L) &&
            isequal(a.P, b.P) 
end

function isequal(a::CMPS, b::CMPS)
    return isequal(a.Q, b.Q) && 
            isequal(a.R, b.R)
end

"""
    Determine if a CMPO is hermitian
"""
ishermitian(o::CMPO) = isequal(o, adjoint(o))


"""
    Perform a unitary transformation in the imaginary-time direction
    (i.e. One physical bond of cMPO and virtual bond of CMPS)
    if U is a square matrix, this is a guage transformation
"""
function project(s::CMPS, u::AbstractMatrix)
    Q = u' * s.Q * u
    if length(size(s.R)) == 2
        R = u' * s.R * u
    else
        R = ein"(ip,pql),qj -> ijl"(u', s.R, u)
    end
    return CMPS(Q, R)
end


function project(o::CMPO, u::AbstractMatrix)
    Q = u' * o.Q * u
    if length(size(o.R)) == 2
        R = u' * o.R * u
        L = u' * o.R * u
        P = u' * o.P * u
    else
        R = ein"(ip,pql),qj -> ijl"(u', o.R, u)
        L = ein"(ip,pql),qj -> ijl"(u', o.L, u)
        P = ein"(ip,pqlm),qj -> ijlm"(u', o.P, u)
    end
    return CMPO(Q,R,L,P)
end


""" multiply a matrix to the left/right of the CMPS/CMPO
    --        --   --            --
    | 1 0 ... 0|   | I + dtau Q   |
    | 0        |   |              |
    | :        |   |       |      |
    | :    W   |   | sqrt(dtau) R |
    | 0        |   |       |      |
    --        --   --            --
"""
function *(Ut::AbstractMatrix, s::CMPS)
    length(size(s.R)) == 2 ? R = Ut[1,1] * s.R :
        R = ein"mn, ijn -> ijm"(Ut, s.R)
    return CMPS(s.Q, R)
end

function *(Ut::AbstractMatrix, o::CMPO)
    if length(size(o.R)) == 2
        R = Ut[1,1] * o.R
        P = Ut[1,1] * o.P
    else
        R = ein"mk, ijk -> ijm"(Ut, o.R)
        P = ein"mk, ijkn -> ijmn"(Ut, o.P)
    end
    return CMPO(o.Q, R, o.L, P)
end

function *(o::CMPO, Ut::AbstractMatrix)
    if length(size(o.L)) == 2
        L = Ut[1,1] * o.L
        P = Ut[1,1] * o.P
    else
        L = ein"ijk, km -> ijm"(o.L, Ut)
        P = ein"ijnk, km -> ijnm"(o.P, Ut) 
    end
    return CMPO(o.Q, o.R, L, P)
end


""" 
    transform the CMPS to the gauge where Q is a diagonalized matrix
"""
function diagQ(s::CMPS)
    _, u = symeigen(s.Q)
    return project(s, u)
end


#end    #module cMPSOperations
