"""
    ChainRules.rrule(::typeof(logtrexp), t::Real, M::AbstractArray, trace_estimator::TraceEstimator{<:typeof{simple_ftlm}})

Rrule for ``ln[Tr(e^{tM})]`` function where `typeof(M) <: CMPSMatrix`.
`∂y_∂M` is also of type `CMPSMatrix`. `trace_estimator = simple_ftlm`.
"""
function ChainRules.rrule(::typeof(logtrexp), t::Real, M::CMPSMatrix, 
                          trace_estimator::TraceEstimator{<:typeof(simple_ftlm), To}) where {To}
    @unpack options = trace_estimator
    @unpack distr, Nr, Nk, processor = options
    @unpack ψl, ψr = M
    Ns = size(M, 1)
    Nk = min(Nk, Ns)

    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)
    solver = solver_function(processor)
    
    sign(t) == 1 ? which = :LR : which=:SR
    processor == CPU ? x0 = rand(T,Ns) : x0 = CUDA.rand(T, Ns)
    e0, _, _ = eigsolve(M, x0, 1, which, ishermitian = true, tol=1.e-12)
    e1 = e0[1]
    expr_Λ = e -> exp(t*(e-e1))
    expr_∂y_∂t = e -> e * exp(t*(e-e1))

    if processor == CPU
        ∂y_∂Ql, ∂y_∂Qr = zeros(T, size(ψl.Q)), zeros(T, size(ψr.Q))
        ∂y_∂Rl, ∂y_∂Rr = zeros(T, size(ψl.R)), zeros(T, size(ψr.R))
    else
        ∂y_∂Ql, ∂y_∂Qr = CUDA.zeros(T, size(ψl.Q)), CUDA.zeros(T, size(ψr.Q))
        ∂y_∂Rl, ∂y_∂Rr = CUDA.zeros(T, size(ψl.R)), CUDA.zeros(T, size(ψr.R))
    end
    y, ∂y_∂t = 0. , 0.
    for r = 1: Nr
        init_vector = random_unit_vector(Ns, distr) |> solver
        @unpack weight, values, vectors = fullortho_lanczos(M; init_vector, Nk) |> eigensolver
        
        y += _sum_expr_w1_w2(expr_Λ, values, weight, conj(weight))
        ∂y_∂t += _sum_expr_w1_w2(expr_∂y_∂t, values, weight, conj(weight))

        Λ = map(expr_Λ, values)
        vectors = reshape(vectors, χr, χl, Nk)
        init_vector = reshape(init_vector, χr, χl)

        ∂y_∂Ql_temp = -t * ein"n,n,kbn,ka,kk -> ab"(Λ, weight, conj(vectors), init_vector)
        ∂y_∂Qr_temp = -t * ein"n,n,bkn,ak,kk -> ab"(Λ, weight, conj(vectors), init_vector)
        ∂y_∂Ql = map(+, ∂y_∂Ql, ∂y_∂Ql_temp)
        ∂y_∂Qr = map(+, ∂y_∂Qr, ∂y_∂Qr_temp)

        if typeof(ψr.R) <: AbstractMatrix
            ∂y_∂Rl_temp = -t * ein"n,n,kl,lbn,ka -> ab"(Λ, weight, ψr.R, conj(vectors), init_vector)
            ∂y_∂Rr_temp = -t * ein"n,n,kl,bln,ak -> ab"(Λ, weight, ψl.R, conj(vectors), init_vector)
        else
            ∂y_∂Rl_temp = -t * ein"n,n,klm,lbn,ka -> abm"(Λ, weight, ψr.R, conj(vectors), init_vector)
            ∂y_∂Rr_temp = -t * ein"n,n,klm,bln,ak -> abm"(Λ, weight, ψl.R, conj(vectors), init_vector)
        end
        ∂y_∂Rl = map(+, ∂y_∂Rl, ∂y_∂Rl_temp)
        ∂y_∂Rr = map(+, ∂y_∂Rr, ∂y_∂Rr_temp)
    end
    
    factor = Ns/Nr
    y = y * factor
    ∂y_∂t = ∂y_∂t * factor /y
    ∂y_∂Ql = ∂y_∂Ql * factor /y
    ∂y_∂Qr = ∂y_∂Qr * factor /y
    ∂y_∂Rl = ∂y_∂Rl * factor /y
    ∂y_∂Rr = ∂y_∂Rr * factor /y

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