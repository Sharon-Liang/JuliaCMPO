import Base: eltype, size, length, getindex, iterate

@with_kw struct CMPSMatrix{T<:Number}
    ψl::AbstractCMPS{T}
    ψr::AbstractCMPS{T}
end

Matrix(A::CMPSMatrix) = A.ψl * A.ψr

eltype(A::CMPSMatrix{T}) where T = T

size(A::CMPSMatrix) = map(*, size(A.ψl.Q), size(A.ψr.Q))
size(A::CMPSMatrix, n) = size(A)[n]
length(A::CMPSMatrix) = size(A,1) * size(A,2)


function getindex(A::CMPSMatrix{T} ,Id::Vararg{Int, 2}) where T
    Nr = size(A.ψr.Q)[1]
    rl, rr = divrem(Id[1]-1, Nr); rl += 1; rr += 1
    cl, cr = divrem(Id[2]-1, Nr); cl += 1; cr += 1
    
    @unpack ψl, ψr = A
    li = convert(typeof(ψl.Q), Matrix{T}(I,size(ψl.Q)))
    ri = convert(typeof(ψr.Q), Matrix{T}(I,size(ψr.Q)))
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)
    res = map(i->L[rl,cl,i]*R[rr,cr,i], 1:size(L)[3]) |> sum
    return -res
end

function getindex(A::CMPSMatrix,i::Int)
    c, r = divrem(i-1, size(A)[2])
    c += 1; r += 1
    return getindex(A, c, r)
end

function *(A::CMPSMatrix{T}, v::AbstractVector{T}) where T
    @unpack ψl, ψr = A
    V = reshape(v, size(ψr.Q)[2], size(ψl.Q)[2])
    V = convert(typeof(ψl.Q), V)

    li = convert(typeof(ψl.Q), Matrix{T}(I,size(ψl.Q)))
    ri = convert(typeof(ψr.Q), Matrix{T}(I,size(ψr.Q)))
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)

    Kv = ein"(ijm,klm),lj -> ki"(L, R, V)
    return vec(-Kv)
end

function *(v::LinearAlgebra.Adjoint{T, Vector{T}}, A::CMPSMatrix{T}) where T
    @unpack ψl, ψr = A
    V = reshape(v, size(ψr.Q)[2], size(ψl.Q)[2])
    V = convert(typeof(ψl.Q), V)

    li = convert(typeof(ψl.Q), Matrix{T}(I,size(ψl.Q)))
    ri = convert(typeof(ψr.Q), Matrix{T}(I,size(ψr.Q)))
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)

    Kv = ein"ki, (ijm,klm) -> lj"(V, L, R)
    return vec(-Kv) |> transpose
end

iterate(A::CMPSMatrix) = (getindex(A,1), 2)
iterate(A::CMPSMatrix, n) = (getindex(A,n), n+1)

#make CMPSMatrix callable
(A::CMPSMatrix)(x) =  A * x

    