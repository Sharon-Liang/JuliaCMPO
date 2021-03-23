#module Setup
#using LinearAlgebra
"""
    D: Virtual bond dimension (depends on Hamitonian)
    χ: Imaginary-time bond dimension
    d: Physical bond dimension （degrees of freedom on one Lattice site)
"""

"""
cMPS -        -
     | 1 + ϵQ |
     |  √ϵR   |
     -        -
     Q: χ × χ matrix
     R: χ × χ × (D-1) array
"""
"""
cMPO -              -
     | 1 + ϵQ   √ϵL |
     |  √ϵR      P  |
     -              -
     Q: d × d matrix : Onsite terms
     R: d × d × (D-1) array   : NN interaction terms
     L: d × d × (D-1) array   : NN interaction terms
     R: d × d × (D-1) × (D-1) array : long-range interaction terms
"""

struct cmps{T<:AbstractArray}
    Q::T
    R::T
end

struct cmpo{T<:AbstractArray}
    Q::T  # onsite
    R::T  # interaction
    L::T  # interaction
    P::T  # long-range
end

function toarray(ψ::cmps)
    # size(Q) == size(R)
    (r,c) = size(ψ.Q)
    x = zeros(r,c,2)
    x[:,:,1] = ψ.Q
    x[:,:,2] = ψ.R
    return x
end

"""multiplications of cmps and cmpo"""
function *(sl::cmps, sr::cmps)
    li = Matrix(1.0I,size(sl.Q))
    ri = Matrix(1.0I,size(sr.Q))
    K = li ⊗ sr.Q + sl.Q ⊗ ri + sl.R ⊗ sr.R
    return -K
end

function *(o::cmpo, s::cmps)
    oi = Matrix(1.0I,size(o.Q))
    si = Matrix(1.0I,size(s.Q))
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return cmps(Q, R)
end

function *(s::cmps, o::cmpo)
    oi = Matrix(1.0I,size(o.Q))
    si = Matrix(1.0I,size(s.Q))
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return cmps(Q, R)
end

function *(ol::cmpo, or::cmpo)
    li = Matrix(1.0I,size(ol.Q))
    ri = Matrix(1.0I,size(or.Q))
    Q = ol.Q ⊗ ri + li ⊗ or.Q + ol.L ⊗ or.R
    L = li ⊗ or.L + ol.L ⊗ or.P
    R = ol.R ⊗ ri + ol.P ⊗ or.R
    P = ol.P ⊗ or.P
    return cmpo(Q,R,L,P)
end

""" init cmps """
function init_cmps(χ::Integer; hermition = true)
    Q = rand(χ, χ)
    R = rand(χ, χ)
    if hermition
        Q = symmetrize(Q)
        R = symmetrize(R)
    end
    return cmps(Q,R)
end

function init_cmps(χ::Integer, W::cmpo)
    # r = 0 case
    d = size(W.Q)[1];  (q,r) = divrem(log(d,χ), 1)
    ψ = cmps(W.Q, W.R)
    if r == 0
        for i = 1:Integer(q-1)  ψ = W * ψ  end
    else
        error("Not support yet :)")
    end
    return ψ
end

function normalize(s::cmps, β::Real)
    T = -β*(s*s) |> trexp
    λ = (T.max + log(T.res))/β
    λ = λ/2
    eye = Matrix(λ*I,size(s.Q))
    Q = s.Q - eye
    return cmps(Q, s.R)
end

function ovlp(sl::cmps, sr::cmps, β::Real)
    -β*(sl*sr) |> trexp |> value
end

function ovlp(s::cmps, β::Real)
    -β*(s*s) |> trexp |> value
end

"""
#end # module
