using OMEinsum
using Zygote, FiniteDiff, ForwardDiff
using LinearAlgebra, GenericLinearAlgebra
import Base: *
using StatsFuns # logsumexp

struct CMPS
    Q::Array{<:Number}
    R::Array{<:Number}
end

struct CMPO
    Q::Array{<:Number}  # onsite
    R::Array{<:Number}  # interaction, column vector
    L::Array{<:Number}  # interaction, row vector
    P::Array{<:Number}  # long-range
end

"""
utilities
"""
function pauli(symbol::Symbol)
    if symbol==:x return [0. 1.; 1. 0.]
    elseif symbol==:y return [0. -1im; 1im 0.]
    elseif symbol==:z return [1. 0.; 0. -1.]
    elseif symbol==:+ return [0. 1.; 0. 0.]
    elseif symbol==:- return [0. 0.; 1. 0.]
    else
        error("The input should be :x,:y,:z,:+,:-.")
    end
end

function TFIsing(J::Real, Γ::Real; field::Symbol=:N, η::Float64 = 1.e-2)
    if field == :N
        h = zeros(2,2)
    else
        h = η .* pauli(field)
    end
    return cmpo(Γ*pauli(:x)+h, √J*pauli(:z), √J*pauli(:z), zeros(2,2))
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

"""
define products of cmps and cmpo
"""
function ⊗(A::Matrix{<:Number}, B::Matrix{<:Number})
    (r1, c1) = size(A)
    (r2,c2) = size(B)
    return reshape(ein"ij,kl->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::Matrix{<:Number}, B::Array{<:Number,3})
    (r1,c1) = size(A)
    (r2,c2,d) = size(B)
    return reshape(ein"ij,klm->kiljm"(A, B), r1*r2, c1*c2, d)
end

function ⊗(A::Array{<:Number,3}, B::Matrix{<:Number})
    (r1,c1,d) = size(A)
    (r2,c2) = size(B)
    return reshape(ein"ijm,kl->kiljm"(A, B), r1*r2, c1*c2, d)
end

function ⊗(A::Array{<:Number,3}, B::Array{<:Number,3})
    (r1,c1,d1) = size(A)
    (r2,c2,d2) = size(B)
    if d1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijm,klm->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::Array{<:Number,4}, B::Array{<:Number,3})
    (r1,c1,d1,f) = size(A)
    (r2,c2,d2) = size(B)
    if f != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijnm,klm->kiljn"(A, B), r1*r2, c1*c2, d1)    
end

function ⊗(A::Array{<:Number,3}, B::Array{<:Number,4}) 
    (r1,c1,d1) = size(A)
    (r2,c2,d2,f) = size(B)
    if d1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijm,klmn->kiljn"(A, B), r1*r2, c1*c2, f)    
end

function ⊗(A::Array{<:Number,4}, B::Array{<:Number,4})
    (r1,c1,d1,f1) = size(A)
    (r2,c2,d2,f2) = size(B)
    if f1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijpm,klmq->kiljpq"(A, B), r1*r2, c1*c2, d1, f2)    
end

function *(sl::CMPS, sr::CMPS)
    li = Matrix(1.0I,size(sl.Q))
    ri = Matrix(1.0I,size(sr.Q))
    K = li ⊗ sr.Q + sl.Q ⊗ ri + sl.R ⊗ sr.R
    return -K
end

function *(o::CMPO, s::CMPS)
    oi = Matrix(1.0I,size(o.Q))
    si = Matrix(1.0I,size(s.Q))
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return CMPS(Q, R)
end

function *(s::CMPS, o::CMPO)
    oi = Matrix(1.0I,size(o.Q))
    si = Matrix(1.0I,size(s.Q))
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return CMPS(Q, R)
end

"""
conversion of cmps and array
"""
function toarray(ψ::CMPS)
    sq = size(ψ.Q)
    sr = size(ψ.R) #sr: dimension of ψ.R array
    if length(sr) == 2 #sr = 2, ψ.R is a matrix
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

function tocmps(A::Array{<:Number,3})
    d = size(A)[3]
    if d == 2
        return CMPS(A[:,:,1],A[:,:,2])
    else
        return CMPS(A[:,:,1],A[:,:,2:end])
    end
end

function symmetrize(A::Matrix{<:Number})
    (A + A')/2
end

function logtrexp(A::Matrix{<:Number})
    #if ishermitian(A) == false
    if isapprox(A,A') == false
        error("The input matrix should be hermitian")
    end
    A = symmetrize(A) |> Hermitian
    return eigvals(A) |> logsumexp
end

function free_energy(ψ::CMPS, W::CMPO, β::Real)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    H = ψ * ψ |> symmetrize |> Hermitian
    res = logtrexp(-β*K)- logtrexp(-β*H)
    return -1/β * res
end

function free_energy(param::Array{<:Number,3}, W::CMPO, β::Real)
    free_energy(tocmps(param), W, β)
end


"""
test begin: compare gradient
"""
function test2(A::AbstractArray)
    ψ = A |> tocmps  
    K = ψ * w * ψ
    return sum(K)
end

function test3(param::Array{T,3} where T<:Number)
    free_energy(tocmps(param), w, 20)
end

w = TFIsing(1.0,1.0)
arr = init_cmps(2,w) |> toarray``
hessian(test2, arr)

# ChainRules: Jin-Guo Liu



# hessian_reverse function is copied from:
# https://github.com/FluxML/Zygote.jl/blob/master/src/lib/grad.jl
"""
    hessian_reverse(f, x)

This should be equivalent to [`hessian(f, x)`](@ref hessian),
but implemented using reverse over reverse mode, all Zygote.
(This is usually much slower, and more likely to find errors.)
"""
#hessian_reverse(f, x::AbstractArray) = jacobian(x -> gradient(f, x)[1], x)[1]
#hessian_reverse(f, x::Number) = gradient(x -> gradient(f, x)[1], x)[1]

function ⊗(A::AbstractMatrix, B::AbstractMatrix)
    (r1, c1) = size(A)
    (r2,c2) = size(B)
    return reshape(ein"ij,kl->kilj"(A, B), r1*r2, c1*c2)
end


function test4(A::AbstractArray)
    ψ = A |> tocmps  
    K = ψ * w * ψ
    return sum(K)
end



function logtrexp(A::AbstractMatrix)
    A = symmetrize(A) |> Hermitian
    eigvals(A) |> logsumexp
end

function test1(A::AbstractMatrix)
    C = A ⊗ A |> symmetrize |> Hermitian
    e = eigvals(C)
    return sum(e)
end

function test2(A::AbstractMatrix)
    return logtrexp(A ⊗ A)
end

function test3(A::AbstractArray)
    ψ = A |> tocmps
    K = ψ * w * ψ
    return logtrexp(K)
end



χ = 2; β = 20
w = TFIsing(1.0,1.0)
arr = init_cmps(χ,w) |> toarray

f1 = x -> free_energy(x, w, β)


K = ψ0 * w * ψ0
hessian(f1, arr)