include("toolfunctions.jl")

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

function toarray(ψ::cMPS)
    # size(Q) == size(R)
    (r,c) = size(ψ.Q)
    x = zeros(r,c,2)
    x[:,:,1] = ψ.Q
    x[:,:,2] = ψ.R
    return x
end

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
    return tr_exp(prod)
end


function difference(ψ1::cMPS, ψ2::cMPS; β=1)
    res = myinnerprod(ψ1,ψ1,β) + myinnerprod(ψ2,ψ2,β)
    res -= myinnerprod(ψ2,ψ1,β) + myinnerprod(ψ1,ψ2,β)
    return res
end
