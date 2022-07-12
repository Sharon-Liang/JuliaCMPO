@reexport import KrylovKit.eigsolve

eigsolve(M::CMPSMatrix{T}, howmany, which; args...) where {T} =
    eigsolve(M, size(M,1), howmany, which, T, ishermitian = true, args...)

    