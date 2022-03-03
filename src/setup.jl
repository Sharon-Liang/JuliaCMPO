#module setup
"""
bD::Int64 # bD: bond dimension along imaginary time direction
vD::Int64 # vD+1: virtual bond dimension of the horizontal legs
pD::Int64 #pD: bond dimension of physical legs
"""
struct CMPS{T<:Number}
    Q::Matrix{T} # Dimension: bD × bD
    R::Array{T}  # Dimension: bD × bD × vD
end

struct CMPO{T<:Number}
    Q::Matrix{T} # Dimension: pD × pD
    R::Array{T}  # Column : pD × pD × vD
    L::Array{T}  # Row: pD × pD × vD
    P::Array{T}  # long-range interaction: pD × pD × vD × vD
end


"""
    Flattern cMPS: cMPS <-> Array
"""
function toarray(ψ::CMPS)
    sq = size(ψ.Q)
    sr = size(ψ.R)  #sr: dimension of ψ.R array
    if length(sr) == 2  #sr = 2, ψ.R is a matrix
        Q = reshape(ψ.Q, sq[1],sq[2],1)
        R = reshape(ψ.R, sr[1],sr[2],1)
    elseif length(sr) > 2
        Q = reshape(ψ.Q, sq[1],sq[2],1)
        R = ψ.R
    else
        @error "Improper CMPS"
    end
    return cat(Q,R,dims=3)
end

function tovector(ψ::CMPS)
    arr = ψ |> toarray
    dim = size(arr)
    return vec(arr), dim
end

function tocmps(A::AbstractArray{T,3} where T)
    d = size(A)[3]
    if d == 2
        return CMPS(A[:,:,1],A[:,:,2])
    else
        return CMPS(A[:,:,1],A[:,:,2:end])
    end
end

function tocmps(V::Vector{T} where T, dim::Tuple)
    arr = reshape(V, dim)
    return tocmps(arr)
end



"""
    Normalize a cMPS, i.e. ⟨ψ|ψ⟩ = 1
"""
function normalize(s::CMPS, β::Real)
    T = -β*(s*s) |> trexp
    λ = (T.max + log(T.res))/β
    λ = λ/2
    eye = Matrix(λ*I,size(s.Q))
    Q = s.Q - eye
    return CMPS(Q, s.R)
end


"""
    Logarithm of the overlap of two cMPS:
        log_overlap = ln(⟨ψl|ψr⟩)
"""
function log_overlap(sl::CMPS, sr::CMPS, β::Real)
    K = sl * sr
    return logtrexp(-β * K)
end


"""
    transpose of a CMPO
"""
function transpose(o::CMPO)
    Q = o.Q
    R = o.L
    L = o.R
    P = ein"ijkl->ijlk"(o.P)
    return CMPO(Q,R,L,P)
end


"""
    Perform a unitary transformation in the imaginary-time direction
    (i.e. One physical bond of cMPO and virtual bond of cMPS)
    if U is a square matrix, this is a guage transformation
"""
function project(s::CMPS, u::AbstractMatrix)
    Q = u' * s.Q * u
    R = ein"(ip,pql),qj -> ijl"(u', s.R, u)
    return CMPS(Q, R)
end

function project(o::CMPO, u::AbstractMatrix)
    Q = u' * o.Q * u
    R = ein"(ip,pql),qj -> ijl"(u', o.R, u)
    L = ein"(ip,pql),qj -> ijl"(u', o.L, u)
    P = ein"(ip,pqlm),qj -> ijlm"(u', o.R, u)
    return CMPO(Q,R,L,P)
end


""" 
    transform the cMPS to the gauge where Q is a diagonalized matrix
"""
function diagQ(s::CMPS; dim::Integer = 999)
    χ = min(size(s.Q)[1], dim)
    _, u = eigensolver(s.Q)
    return project(s, u[:,end+1-χ:end])
end


#end module setup
