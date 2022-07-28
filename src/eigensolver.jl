@reexport import KrylovKit.eigsolve
@reexport import FiniteTLanczos.eigensolver

"""
    eigensolver(A::CMPSMatrix)

Solve all eigen pairs of a symmetric `A`.
"""
function eigensolver(A::CMPSMatrix)
    res = tomatrix(A) |> symmetrize
    return eigensolver(res)
end
