#module utilities
import Base: kron

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

function delta(x::Real, η::Real)
    num = η
    den = (x^2 + η^2) * π
    return num/den
end

function ⊗(A::AbstractMatrix, B::AbstractMatrix)
    (r1, c1) = size(A)
    (r2,c2) = size(B)
    return reshape(ein"ij,kl->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::AbstractMatrix, B::Array{Float64,3})
    (r1,c1) = size(A)
    (r2,c2,d) = size(B)
    return reshape(ein"ij,klm->kiljm"(A, B), r1*r2, c1*c2, d)
end

function ⊗(A::Array{Float64,3}, B::AbstractMatrix)
    (r1,c1,d) = size(A)
    (r2,c2) = size(B)
    return reshape(ein"ijm,kl->kiljm"(A, B), r1*r2, c1*c2, d)
end

function ⊗(A::Array{Float64,3} , B::Array{Float64,3})
    (r1,c1,d1) = size(A)
    (r2,c2,d2) = size(B)
    if d1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijm,klm->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::Array{Float64,4}, B::Array{Float64,3})
    (r1,c1,d1,f) = size(A)
    (r2,c2,d2) = size(B)
    if f != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijnm,klm->kiljn"(A, B), r1*r2, c1*c2, d1)    
end

function ⊗(A::Array{Float64,3}, B::Array{Float64,4})
    (r1,c1,d1) = size(A)
    (r2,c2,d2,f) = size(B)
    if d1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijm,klmn->kiljn"(A, B), r1*r2, c1*c2, f)    
end

function ⊗(A::Array{Float64,4}, B::Array{Float64,4})
    (r1,c1,d1,f1) = size(A)
    (r2,c2,d2,f2) = size(B)
    if f1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijpm,klmq->kiljpq"(A, B), r1*r2, c1*c2, d1, f2)    
end

function symmetrize(A::AbstractMatrix)
    (A + A')/2
end

struct trexp{T<:Number}
    max::T
    res::T
end

function value(x::trexp)
    exp(x.max) * x.res
end

function trexp(A::AbstractMatrix)
    if ishermitian == false
        error("The input matrix should be hermitian")
    end
    A = symmetrize(A) |> Hermitian
    val= eigvals(A)
    max = maximum(val)
    res = exp.(val .- max) |> sum
    trexp(max, res)
end

function logtrexp(A::AbstractMatrix)
    if ishermitian == false
        error("The input matrix should be hermitian")
    end
    A = symmetrize(A) |> Hermitian
    eigvals(A) |> logsumexp
end

"""function manipulation"""
function grad_func(f::Function, var::AbstractArray)
    function gf(gx::AbstractArray, var::AbstractArray)
        gx[1:end] = gradient(var -> f(var),var)[1][1:end]
    end
    gf
end

#end  # module utilities
