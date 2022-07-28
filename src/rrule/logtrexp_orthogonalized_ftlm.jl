"""
    ChainRules.rrule(::typeof(logtrexp), t::Real, M::AbstractArray, trace_estimator::TraceEstimator{<:typeof{orthogonalized_ftlm}})

Rrule for ``ln[Tr(e^{tM})]`` function where `typeof(M) <: CMPSMatrix`.
`∂y_∂M` is also of type `CMPSMatrix`. `trace_estimator = orthogonalized_ftlm`.
"""
function ChainRules.rrule(::typeof(logtrexp), t::Real, M::CMPSMatrix, 
                          trace_estimator::TraceEstimator{<:typeof(orthogonalized_ftlm), To}) where {To}
    @unpack options = trace_estimator
    @unpack distr, Nr, Nk, Ne, processor = options
    @unpack ψl, ψr = M
    solver = solver_function(processor)

    Ns = size(M, 1)
    Ne = min(Ne, Ns)
    Nk = min(Nk, Ns - Ne) 

    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)

    sign(t) == 1 ? which = :LR : which=:SR
    processor == CPU ? x0 = rand(T, Ns) : x0 = CUDA.rand(T, Ns)
    e0, _, _ =eigsolve(M, x0, 1, which, ishermitian = true, tol=1.e-12)
    e1 = e0[1]
    expr_Λ = e -> exp(t*(e-e1))
    expr_∂y_∂t = e -> e * exp(t*(e-e1))

    #exact eigen pair part
    krylovdim = max(30, Ne+10)
    vals, vecs, _ = eigsolve(M, x0, Ne, which, ishermitian = true, krylovdim=krylovdim)
    eigen_vals = solver(vals[1:Ne])
    eigen_vecs = hcat(vecs[1:Ne]...)

    y1 = reduce(+, map(expr_Λ, eigen_vals))
    ∂y_∂t1 = reduce(+, map(expr_Λ, eigen_vals))
    
    Λ1 = map(expr_Λ, eigen_vals)
    eigen_vecs2 = reshape(eigen_vecs, χr, χl, Ne)

    ∂y_∂Ql1 = -t * ein"n,kbn,kan -> ab"(Λ1, conj(eigen_vecs2), eigen_vecs2)
    ∂y_∂Qr1 = -t * ein"n,bkn,akn -> ab"(Λ1, conj(eigen_vecs2), eigen_vecs2)
    if typeof(ψr.R) <: AbstractMatrix
        ∂y_∂Rl1 = -t * ein"n,kl,lbn,kan,kl -> ab"(Λ1, ψr.R, conj(eigen_vecs2), eigen_vecs2)
        ∂y_∂Rr1 = -t * ein"n,kl,bln,akn,kl -> ab"(Λ1, ψl.R, conj(eigen_vecs2), eigen_vecs2)
    else
        ∂y_∂Rl1 = -t * ein"n,klm,lbn,kan,klm -> abm"(Λ1, ψr.R, conj(eigen_vecs2), eigen_vecs2)
        ∂y_∂Rr1 = -t * ein"n,klm,bln,akn,klm -> abm"(Λ1, ψl.R, conj(eigen_vecs2), eigen_vecs2)
    end

    #lanczos part
    y2, ∂y_∂t2 = 0. , 0.
    if processor == CPU
        ∂y_∂Ql2, ∂y_∂Qr2 = zeros(T, size(ψl.Q)), zeros(T, size(ψr.Q))
        ∂y_∂Rl2, ∂y_∂Rr2 = zeros(T, size(ψl.R)), zeros(T, size(ψr.R))
    else
        ∂y_∂Ql2, ∂y_∂Qr2 = CUDA.zeros(T, size(ψl.Q)), CUDA.zeros(T, size(ψr.Q))
        ∂y_∂Rl2, ∂y_∂Rr2 = CUDA.zeros(T, size(ψl.R)), CUDA.zeros(T, size(ψr.R))
    end
    if Nk > 0
        for r = 1: Nr
            init_vector = random_unit_vector(Ns, distr) |> solver
            @unpack weight, values, vectors = fullortho_lanczos(M; init_vector, Nk) |> eigensolver

            y2 += _sum_expr_w1_w2(expr_Λ, values, weight, conj(weight))
            ∂y_∂t2 += _sum_expr_w1_w2(expr_∂y_∂t, values, weight, conj(weight))

            Λ2 = map(expr_Λ, values)
            vectors = reshape(vectors, χr, χl, Nk)
            init_vector = reshape(init_vector, χr, χl)

            ∂y_∂Ql2_temp = -t * ein"n,n,kbn,ka -> ab"(Λ2, weight, conj(vectors), init_vector)
            ∂y_∂Qr2_temp = -t * ein"n,n,bkn,ak -> ab"(Λ2, weight, conj(vectors), init_vector)
            ∂y_∂Ql2 = map(+, ∂y_∂Ql2, ∂y_∂Ql2_temp)
            ∂y_∂Qr2 = map(+, ∂y_∂Qr2, ∂y_∂Qr2_temp)

            if typeof(ψr.R) <: AbstractMatrix
                ∂y_∂Rl2_temp = -t * ein"n,n,kl,lbn,ka -> ab"(Λ2, weight, ψr.R, conj(vectors), init_vector)
                ∂y_∂Rr2_temp = -t * ein"n,n,kl,bln,ak -> ab"(Λ2, weight, ψl.R, conj(vectors), init_vector)
            else
                ∂y_∂Rl2_temp = -t * ein"n,n,klm,lbn,ka -> abm"(Λ2, weight, ψr.R, conj(vectors), init_vector)
                ∂y_∂Rr2_temp = -t * ein"n,n,klm,bln,ak -> abm"(Λ2, weight, ψl.R, conj(vectors), init_vector)
            end
            ∂y_∂Rl2 = map(+, ∂y_∂Rl2, ∂y_∂Rl2_temp)
            ∂y_∂Rr2 = map(+, ∂y_∂Rr2, ∂y_∂Rr2_temp)
        end
    end
    
    factor = (Ns-Ne)/Nr
    y = y1 + y2*factor
    ∂y_∂t = (∂y_∂t1 + ∂y_∂t2 * factor) /y
    ∂y_∂Ql = (∂y_∂Ql1 + ∂y_∂Ql2 * factor) /y
    ∂y_∂Qr = (∂y_∂Qr1 + ∂y_∂Qr2 * factor) /y
    ∂y_∂Rl = (∂y_∂Rl1 + ∂y_∂Rl2 * factor) /y
    ∂y_∂Rr = (∂y_∂Rr1 + ∂y_∂Rr2 * factor) /y

    
    function logtrexp_pullback(ȳ)
        t̄ = ȳ * ∂y_∂t

        Q̄l = ȳ * ∂y_∂Ql
        R̄l = ȳ * ∂y_∂Rl
        Q̄r = ȳ * ∂y_∂Qr
        R̄r = ȳ * ∂y_∂Rr
        ψ̄l = CMPS_generate(Q̄l, R̄l)
        ψ̄r = CMPS_generate(Q̄r, R̄r)
        M̄ = CMPSMatrix(ψ̄l, ψ̄r)
        return ChainRules.NoTangent(), t̄, M̄, ChainRules.NoTangent()
    end
    return y, logtrexp_pullback
end