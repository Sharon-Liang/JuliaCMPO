#module utilities
import Base: kron

"""
    Pauli matrices
"""
function pauli(T::DataType, symbol::Symbol)
    One = one(T)
    Zero = zero(T)
    if symbol==:x return [Zero One; One Zero]
    elseif symbol==:y return [Zero -One*im; One*im Zero]
    elseif symbol==:iy return [Zero One; -One Zero]
    elseif symbol==:z return [One Zero; Zero -One]
    elseif symbol==:+ return [Zero One; Zero Zero]
    elseif symbol==:- return [Zero Zero; One Zero]
    else
        error("The input should be :x,:y,:z,:+,:-, :iy.")
    end
end

function pauli(symbol::Symbol)
    return pauli(Float64, symbol)
end


"""
Gaussian approximate delta function
"""
function delta(x::Real, η::Real)
    num = η
    den = (x^2 + η^2) * π
    return num/den
end


"""
    Masubara frequency
"""
function Masubara_freq(n::Int64, β::Real; type::Symbol=:b)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end


"""
    Products
"""
function ⊗(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    (r1, c1) = size(A)
    (r2, c2) = size(B)
    return reshape(ein"ij,kl->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::AbstractMatrix{T}, B::AbstractArray{T,3}) where T
    (r1, c1) = size(A)
    (r2, c2, d) = size(B)
    return reshape(ein"ij,klm->kiljm"(A, B), r1*r2, c1*c2, d)
end

function ⊗(A::AbstractArray{T,3}, B::AbstractMatrix{T}) where T
    (r1, c1, d) = size(A)
    (r2, c2) = size(B)
    return reshape(ein"ijm,kl->kiljm"(A, B), r1*r2, c1*c2, d)
end

function ⊗(A::AbstractArray{T,3}, B::AbstractArray{T,3}) where T
    (r1, c1, d1) = size(A)
    (r2, c2, d2) = size(B)
    if d1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijm,klm->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::AbstractArray{T,4}, B::AbstractArray{T,3}) where T
    (r1, c1, d1, f) = size(A)
    (r2, c2, d2) = size(B)
    if f != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijnm,klm->kiljn"(A, B), r1*r2, c1*c2, d1)    
end

function ⊗(A::AbstractArray{T,3}, B::AbstractArray{T,4})  where T
    (r1, c1, d1) = size(A)
    (r2, c2, d2, f) = size(B)
    if d1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijm,klmn->kiljn"(A, B), r1*r2, c1*c2, f)    
end

function ⊗(A::AbstractArray{T,4}, B::AbstractArray{T,4}) where T
    (r1,c1,d1,f1) = size(A)
    (r2,c2,d2,f2) = size(B)
    if f1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijpm,klmq->kiljpq"(A, B), r1*r2, c1*c2, d1, f2)    
end


"""
    Symmetrize of matrices
"""
function symmetrize(A::AbstractMatrix)
    (A + A')/2
end


"""
    Tr[exp(A)] function
"""
struct TrExp{T<:Number}
    max::T
    res::T
end

function value(x::TrExp)
    exp(x.max) * x.res
end

function trexp(A::AbstractMatrix)
    #if ishermitian(A) == false
    if isapprox(A,A') == false
        error("The input matrix should be hermitian")
    end
    A = symmetrize(A) |> Hermitian
    val = eigvals(A)
    max = maximum(val)
    res = exp.(val .- max) |> sum
    return TrExp(max, res)
end


"""
    ln(Tr[exp(A)]) function
"""
function logtrexp(A::AbstractMatrix)
    #if isapprox(A,A') == false
    #    error("The input matrix should be hermitian")
    #end
    A = symmetrize(A) |> Hermitian
    return eigvals(A) |> logsumexp
end



"""function manipulation"""
function gradient_function(f::Function)
    function func(val::Array{T}, var::Array{T}) where T<:Number
        val[1:end] = gradient(f,var)[1][1:end]
    end
    return func
end

function hessian_function(f::Function)
    function func(val::Matrix{T}, var::Vector{T}) where T<:Number
        val[1:end] = hessian(f,var)[1:end]
    end
    return func
end

#end  # module utilities
