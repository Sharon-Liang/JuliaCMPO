@reexport import KrylovKit.eigsolve
import LinearAlgebra.Eigen

function eigensolver(M::CMPSMatrix{Ts,T,S,U}; 
                     which::Symbol=:SR, 
                     ishermitian = true, 
                     kwargs...) where {Ts,T,S,U} 
    N = size(M,1)
    Ts <: CMPS ? x0 = rand(T,N) : x0 = CUDA.rand(T, N)
    vals, vecs, _ = eigsolve(M, x0, N, which, ishermitian = ishermitian, krylovdim = N, kwargs...)
    vals = convert(typeof(x0), vals)
    vecs = hcat(vecs[1:N]...)
    return Eigen(vals, vecs)
end

#function eigsolve(M::CMPSMatrix{Ts,T,S,U}, x0, howmany, which; args...) where {Ts,T,S,U}
#    Ts <: CMPS ? x0 = rand(size(M,1)) : x0 = CUDA.rand(size(M,1))
#    return eigsolve(M, x0, howmany, which, T, ishermitian = true, args...)
#end

    