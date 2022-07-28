import Base: eltype, size, length, getindex, iterate
"""
    CMPSMatrix

Memory saving structure for matrix ``⟨ψl|ψr⟩``.
"""
@with_kw struct CMPSMatrix{Ts<:AbstractCMPS, T, S, U}
    ψl::Ts
    ψr::Ts
    CMPSMatrix{Ts,T,S,U}(ψl::AbstractCMPS{T,S,U}, 
        ψr::AbstractCMPS{T,S,U}) where {Ts,T,S,U} = new(ψl, ψr)
end
CMPSMatrix(ψl::AbstractCMPS{T,S,U}, ψr::AbstractCMPS{T,S,U}) where {T,S,U} = 
    CMPSMatrix{typeof(ψl),T,S,U}(ψl, ψr)

"""
    eltype(A::CMPSMatrix)

Element type of a CMPSMatrix `A`.
""" 
eltype(::CMPSMatrix{Ts,T,S,U}) where {Ts,T,S,U} = T


"""
    size(A::CMPSMatrix)
    size(A::CMPSMatrix, n)
""" 
size(A::CMPSMatrix) = map(*, size(A.ψl.Q), size(A.ψr.Q))
size(A::CMPSMatrix, n) = size(A)[n]


"""
    length(A::CMPSMatrix)
""" 
length(A::CMPSMatrix) = size(A,1) * size(A,2)


"""
    getindex(A::CMPSMatrix, Id::Vararg{Int, 2})
    getindex(A::CMPSMatrix, i::Int)
""" 
function getindex(A::CMPSMatrix ,Id::Vararg{Int, 2})
    Nr = size(A.ψr.Q, 1)
    rl, rr = divrem(Id[1]-1, Nr); rl += 1; rr += 1
    cl, cr = divrem(Id[2]-1, Nr); cl += 1; cr += 1
    
    @unpack ψl, ψr = A
    li, ri = oneunit(ψl.Q), oneunit(ψr.Q)
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)
    return reduce(+, L[rl,cl,:] .* R[rr,cr,:])
end

function getindex(A::CMPSMatrix, i::Int)
    c, r = divrem(i-1, size(A,2))
    c += 1; r += 1
    return getindex(A, c, r)
end


"""
    *(A::CMPSMatrix, v::AbstractVector)

Compute the matrix-vector product of a CMPSMatrix and a vector.
""" 
function *(A::CMPSMatrix, v::AbstractVector)
    @unpack ψl, ψr = A
    V = reshape(v, size(ψr.Q,2), size(ψl.Q,2))
    V = convert(typeof(ψl.Q), V)

    li, ri = oneunit(ψl.Q), oneunit(ψr.Q)
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)

    Kv = ein"(ijm,klm),lj -> ki"(L, R, V)
    return vec(-Kv)
end

function *(v::LinearAlgebra.Adjoint{T, Vector{T}}, A::CMPSMatrix) where {T}
    @unpack ψl, ψr = A
    V = reshape(v, size(ψr.Q,2), size(ψl.Q,2))
    V = convert(typeof(ψl.Q), V)

    li, ri = oneunit(ψl.Q), oneunit(ψr.Q)
    L = cat(li, ψl.Q, ψl.R, dims = 3)
    R = cat(ψr.Q, ri, ψr.R, dims = 3)

    Kv = ein"ki, (ijm,klm) -> lj"(V, L, R)
    return vec(-Kv) |> transpose
end


"""
    iterate(A::CMPSMatrix)
    iterate(A::CMPSMatrix, n::Int)
""" 
iterate(A::CMPSMatrix) = (getindex(A,1), 2)
iterate(A::CMPSMatrix, n) = (getindex(A,n), n+1)

#make CMPSMatrix callable
(A::CMPSMatrix)(x) =  A * x

    
"""
    CTensor(x::CMPSMatrix)

Convert tensors in `CMPSMatrix` to `CMPS`.
"""
CTensor(x::CMPSMatrix) = CMPSMatrix(CTensor(x.ψl), CTensor(x.ψr))

"""
    CuCTensor(x::CMPSMatrix)

Convert tensors in `CMPSMatrix` to `CuCMPS`.
"""
CuCTensor(x::CMPSMatrix) = CMPSMatrix(CuCTensor(x.ψl), CuCTensor(x.ψr))

