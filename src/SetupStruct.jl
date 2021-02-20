module SetupStruct

using LinearAlgebra
using Reexport: @reexport

include("ToolFunctions.jl")
@reexport using .ToolFunctions

export cmps, cmpo
export toarray, myprod, myinnerprod
export difference

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
struct cmps{T<:AbstractArray}
    Q::T
    R::T
end

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

function myprod(O::cmpo, S::cmps)
    Oi = Matrix(1.0I,size(O.Q))
    Si = Matrix(1.0I,size(S.Q))
    Q = kron(Oi , S.Q) + kron(O.Q , Si) + kron(O.L , S.R)
    R = kron(O.R , Si) + kron(O.P , S.R)
    return cmps(Q, R)
end

function myinnerprod(sl::cmps, sr::cmps, β::Real)
    li = Matrix(1.0I,size(sl.Q))
    ri = Matrix(1.0I,size(sr.Q))
    prod -= kron(li , sr.Q) + kron(sl.Q , ri) + kron(sl.R, sr.R)
    return tr_exp(prod)
end


function difference(ψ1::cmps, ψ2::cmps; β=1)
    res = myinnerprod(ψ1,ψ1,β) + myinnerprod(ψ2,ψ2,β)
    res -= myinnerprod(ψ2,ψ1,β) + myinnerprod(ψ1,ψ2,β)
    return res
end

"""
struct option 2:
    struct cMPS(Ts<:AbstractArray, Td<:Integer)
        χ::Td
        D::Td
        Q::Ts
        R::Ts
        cmps::Ts
        fuction make_cmps(Q,R)
            cmps = zeros(χ, χ, D)
            for i = 1:χ, j = 1:χ
                cmps[i,j,1]= Q[i,j]
                for k = 1:D-1  cmps[i,j,k+1] = R[i,j,k] end
            end
        end
    end
pro: 1. natural input for optim.optimizer
     2. products: Einsum
        (mind the order: ij,kl -> kilj |> matrix |> tr_exp)
con: complicated cMPO setup (avoid Array{Array,N}, （D-1）× (D-1) or D^2 × 1 ?)
"""
end  # module SetupStruct
