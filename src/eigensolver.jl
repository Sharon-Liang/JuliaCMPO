import FiniteTLanczos.eigensolver
import KrylovKit.eigsolve

""" eigensolver for CuMatrix
    For CUDA dense matrix eigen solver `CUSOLVER.syevd!` and `CUSOLVER.heevd!`:
        'N'/'V': return eigenvalues/both eigenvalues and eigenvectors 
        'U'/'L': Upper/Lower triangle of `M` is stored.
"""
for elty in (:Float32, :Float64)
    @eval begin
        eigensolver(M::CuMatrix{$elty}) = 
            CUSOLVER.syevd!('V', 'U', M)
    end
end

for elty in (:ComplexF32, :ComplexF64)
    @eval begin
        eigensolver(M::CuArray{$elty}) = 
            CUSOLVER.heevd!('V', 'U', M)
    end
end

eigsolve(M::CMPSMatrix{T}, howmany, which; args...) where {T} =
    eigsolve(M, size(M,1), howmany, which, T, ishermitian = true, args...)

    