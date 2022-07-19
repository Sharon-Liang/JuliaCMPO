@reexport import KrylovKit.eigsolve
import LinearAlgebra.Eigen

function eigensolver(M::CMPSMatrix{Ts,T,S,U}) where {Ts,T,S,U} 
    res = Matrix(M) |> symmetrize
    return eigensolver(res)
end

#end