@reexport import KrylovKit.eigsolve

function eigsolve(M::CMPSMatrix{Ts,T,S,U}, howmany, which; args...) where {Ts,T,S,U}
    Ts <: CMPS ? x0 = rand(size(M,1)) : x0 = CUDA.rand(size(M,1))
    return eigsolve(M, x0, howmany, which, T, ishermitian = true, args...)
end

    