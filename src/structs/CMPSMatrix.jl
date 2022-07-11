import Base: eltype, size, length, getindex, iterate
"""
    `CMPSMatrix(ψl,ψr)`: memory saving struct of matrix `⟨ψl|ψr⟩`
"""
@with_kw struct CMPSMatrix{Ts <: AbstractCMPS, T, S, U}
    ψl::Ts
    ψr::Ts
    CMPSMatrix{Ts,T,S,U}(ψl::AbstractCMPS{T,S,U}, 
        ψr::AbstractCMPS{T,S,U}) where {Ts,T,S,U} = new(ψl, ψr)
end
CMPSMatrix(ψl::AbstractCMPS{T,S,U}, ψr::AbstractCMPS{T,S,U}) where {T,S,U} = 
    CMPSMatrix{typeof(ψl),T,S,U}(ψl, ψr)

"""
    `eltype` of CMPSMatrix
""" 
eltype(A::CMPSMatrix{Ts,T,S,U}) where {Ts,T,S,U} = T

"""
    `size` of CMPSMatrix
""" 
size(A::CMPSMatrix) = map(*, size(A.ψl.Q), size(A.ψr.Q))
size(A::CMPSMatrix, n) = size(A)[n]


"""
    `length` of CMPSMatrix
""" 
length(A::CMPSMatrix) = size(A,1) * size(A,2)


"""
    `getindex` of CMPSMatrix
""" 
function getindex(A::CMPSMatrix{Ts,T,S,U} ,Id::Vararg{Int, 2}) where {Ts,T,S,U}
    Nr = size(A.ψr.Q)[1]
    rl, rr = divrem(Id[1]-1, Nr); rl += 1; rr += 1
    cl, cr = divrem(Id[2]-1, Nr); cl += 1; cr += 1
    
    @unpack ψl, ψr = A
    li = convert(S, Matrix{T}(I,size(ψl.Q)))
    ri = convert(S, Matrix{T}(I,size(ψr.Q)))
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)
    res = map(i->L[rl,cl,i]*R[rr,cr,i], 1:size(L)[3]) |> sum
    return -res
end

function getindex(A::CMPSMatrix,i::Int)
    c, r = divrem(i-1, size(A,2))
    c += 1; r += 1
    return getindex(A, c, r)
end


"""
    Matrix-vector product of a CMPSMatrix and a vector
""" 
function *(A::CMPSMatrix{Ts,T,S,U}, v::AbstractVector{T}) where {Ts,T,S,U}
    @unpack ψl, ψr = A
    V = reshape(v, size(ψr.Q,2), size(ψl.Q,2))
    V = convert(S, V)

    li = convert(S, Matrix{T}(I,size(ψl.Q)))
    ri = convert(S, Matrix{T}(I,size(ψr.Q)))
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)

    Kv = ein"(ijm,klm),lj -> ki"(L, R, V)
    return vec(-Kv)
end

function *(v::LinearAlgebra.Adjoint{T, Vector{T}}, A::CMPSMatrix{Ts,T,S,U}) where {Ts,T,S,U}
    @unpack ψl, ψr = A
    V = reshape(v, size(ψr.Q,2), size(ψl.Q,2))
    V = convert(S, V)

    li = convert(S, Matrix{T}(I,size(ψl.Q)))
    ri = convert(S, Matrix{T}(I,size(ψr.Q)))
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)

    Kv = ein"ki, (ijm,klm) -> lj"(V, L, R)
    return vec(-Kv) |> transpose
end


"""
    `iterate` of CMPSMatrix
""" 
iterate(A::CMPSMatrix) = (getindex(A,1), 2)
iterate(A::CMPSMatrix, n) = (getindex(A,n), n+1)

#make CMPSMatrix callable
(A::CMPSMatrix)(x) =  A * x

    
"""
    convert tensors in CMPSMatrix to CTensor/CuCTensor
"""
CTensor(x::T) where T<:CMPSMatrix = CMPSMatrix(ψl=CTensor(x.ψl), ψr=CTensor(x.ψr))
CuCTensor(x::T) where T<:CMPSMatrix = CMPSMatrix(ψl=CuCTensor(x.ψl), ψr=CuCTensor(x.ψr))