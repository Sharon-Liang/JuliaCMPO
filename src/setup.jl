#module setup
struct CMPS{T<:Number}
    Q::Matrix{T}
    R::Array{T}
end

struct CMPO{T<:Number}
    Q::Matrix{T}  #onsite
    R::Array{T}  #interaction, column vector
    L::Array{T}  #interaction, row vector
    P::Array{T}  #long-range
end

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



"""multiplications of cmps and cmpo"""
function *(sl::CMPS{T}, sr::CMPS{T}) where T
    li = Matrix{T}(I,size(sl.Q))
    ri = Matrix{T}(I,size(sr.Q))
    K = li ⊗ sr.Q + sl.Q ⊗ ri + sl.R ⊗ sr.R
    return -K
end

function *(o::CMPO{T}, s::CMPS{T}) where T
    oi = Matrix(1.0I,size(o.Q))
    si = Matrix(1.0I,size(s.Q))
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return CMPS(Q, R)
end

function *(s::CMPS{T}, o::CMPO{T}) where T
    oi = Matrix{T}(I,size(o.Q))
    si = Matrix{T}(I,size(s.Q))
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return CMPS(Q, R)
end

function *(ol::CMPO{T}, or::CMPO{T}) where T
    li = Matrix{T}(I,size(ol.Q))
    ri = Matrix{T}(I,size(or.Q))
    Q = ol.Q ⊗ ri + li ⊗ or.Q + ol.L ⊗ or.R
    L = li ⊗ or.L + ol.L ⊗ or.P
    R = ol.R ⊗ ri + ol.P ⊗ or.R
    P = ol.P ⊗ or.P
    return CMPO(Q,R,L,P)
end


""" init cmps """
function init_cmps(χ::Int64; D::Int64 = 1, hermition = true, dtype = Float64)
    Q = rand(dtype, χ, χ)
    if D == 1
        R = rand(dtype, χ, χ)
    else
        R = rand(dtype, χ, χ, D)
    end

    if hermition
        Q = symmetrize(Q)
        for d = 1:D
            R[:,:,d] = symmetrize(R[:,:,d])
        end
    end
    return CMPS(Q,R)
end

function init_cmps(χ::Int64, W::CMPO)
    # r = 0 case
    d = size(W.Q)[1];  (q,r) = divrem(log(d,χ), 1)
    ψ = CMPS(W.Q, W.R)
    if r == 0
        for i = 1:Integer(q-1)  ψ = W * ψ  end
    else
        error("Not support yet :)")
    end
    return ψ
end

function normalize(s::CMPS, β::Real)
    T = -β*(s*s) |> trexp
    λ = (T.max + log(T.res))/β
    λ = λ/2
    eye = Matrix(λ*I,size(s.Q))
    Q = s.Q - eye
    return CMPS(Q, s.R)
end

function ovlp(sl::CMPS, sr::CMPS, β::Real)
    -β*(sl*sr) |> trexp |> value
end

function ovlp(s::CMPS, β::Real)
    -β*(s*s) |> trexp |> value
end

#end module setup
