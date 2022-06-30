"""
`ln(Tr[exp(A)])` function for a hermitian matrix `A`
"""
function logtrexp end
    
function logtrexp(A::AbstractMatrix)
    A = symmetrize(A)
    val, _ = eigensolver(A)
    return logsumexp(val)
end

function logtrexp(M::CMPSMatrix{T}, method::Function) where T
    e0, _, _ = eigsolve(M, size(M,1), 1, :LR, ishermitian = true)
    e0 = e0[1]
    expr = e -> exp(e .- e0)
    res = method(M, expr)[1]
    return log(res) + e0
end