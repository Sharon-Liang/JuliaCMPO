@reexport import FiniteTLanczos.eigensolver
@reexport import KrylovKit.eigsolve

for elty in (:Float32, :Float64)
    """ 
    eigensolver for CuMatrix
    For CUDA dense matrix eigen solver `CUSOLVER.syevd!` and `CUSOLVER.heevd!`:
        'N'/'V': return eigenvalues/both eigenvalues and eigenvectors 
        'U'/'L': Upper/Lower triangle of `M` is stored.
    """
    @eval eigensolver(M::CuMatrix{$elty}) = CUSOLVER.syevd!('V', 'U', M)
end

for elty in (:ComplexF32, :ComplexF64)
    """ 
    eigensolver for CuMatrix
    For CUDA dense matrix eigen solver `CUSOLVER.syevd!` and `CUSOLVER.heevd!`:
        'N'/'V': return eigenvalues/both eigenvalues and eigenvectors 
        'U'/'L': Upper/Lower triangle of `M` is stored.
    """
    @eval eigensolver(M::CuArray{$elty}) = CUSOLVER.heevd!('V', 'U', M)
end

eigsolve(M::CMPSMatrix{T}, howmany, which; args...) where {T} =
    eigsolve(M, size(M,1), howmany, which, T, ishermitian = true, args...)

    