using LinearAlgebra
using Zygote
using Optim
using Random
Random.seed!()
"""Tool Functions"""

"""Structs"""
struct cMPS{T<:AbstractArray}
    Q::T
    R::T
end

struct cMPO{T<:AbstractArray}
    Q::T  # onsite
    R::T  # interaction
    L::T  # interaction
    P::T  # long-range
end

"""Functions"""
function myprod(O::cMPO, S::cMPS)
    Oi = Matrix(1.0I,size(O.Q))
    Si = Matrix(1.0I,size(S.Q))
    Q = kron(Oi , S.Q) + kron(O.Q , Si) + kron(O.L , S.R)
    R = kron(O.R , Si) + kron(O.P , S.R)
    return cMPS(Q, R)
end

function myinnerprod(sl::cMPS, sr::cMPS, β::Real)
    li = Matrix(1.0I,size(sl.Q))
    ri = Matrix(1.0I,size(sr.Q))
    prod = kron(li , sr.Q) + kron(sl.Q , ri) + kron(sl.R, sr.R)
    vals = eigvals(prod)
    res = 0.
    for i=1:length(vals)
        res += exp(-β*vals[i])
    end
    return res
end

function F(ψ::cMPS, W::cMPO, β::Real)
    Hψ = myprod(W,ψ)
    res = log(myinnerprod(ψ, Hψ ,β))- log(myinnerprod(ψ,ψ,β))
    return -1/β * res
end

"""Setups"""
J = 1.0; Γ = 1.0
X = [0. 1.; 1. 0.]
Z = [1. 0.; 0. -1.]
β = 1.

W = cMPO(Γ*X, √J*Z, √J*Z, zeros(2,2))
Q = rand(2,2)
R = rand(2,2);
ψ = cMPS(Q,R)

# test
gradient(ψ -> F(ψ, W, β), ψ)
