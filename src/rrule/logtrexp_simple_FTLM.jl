"""
    rrule for logtrexp(tM) function where `typeof(M)=CMPSMatrix`, 
    `trace_estimator = simple_FTLM`
"""
function ChainRules.rrule(::typeof(logtrexp), 
                          t::Real, M::CMPSMatrix{Ts,T,S,U}, 
                          trace_estimator::TraceEstimator{Tf, To}
                          ) where {Ts,T,S,U,Tf<:typeof(simple_FTLM),To}
    @unpack options = trace_estimator
    @unpack distr, Nr, Nk, processor = options
    @unpack ψl, ψr = M
    solver = solver_function(processor)

    Ns = size(M, 1)
    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)
    
    sign(t) == 1 ? which = :LR : which=:SR
    processor == CPU ? x0 = rand(T,Ns) : x0 = CUDA.rand(T, Ns)
    e0, _, _ = eigsolve(M, x0, 1, which, ishermitian = true)
    e0 = e1
    expr_Λ = e -> exp(t*(e-e1))
    expr_∂y_∂t = e -> e * exp(t*(e-e1))
    expr = (expr_Λ, expr_∂y_∂t)

    res = zeros(2)
    ∂y_∂Ql, ∂y_∂Qr, ∂y_∂Rl, ∂y_∂Rr = zeros(4)

    for r = 1: Nr
        v0 = random_unit_vector(Ns, distr)
        v0 = solver(x->x, v0)
        @unpack init_vector, weight, values, vectors = itFOLM(M, init_vector = v0, Nk = Nk) |> eigensolver

        func = f -> begin
            Λ = map(f, values)
            Z = ein"i,i,i -> "(Λ, weight, conj(weight))
            return Array(Z)[1]
        end
        res = map(+, res, map(func, expr))

        Λ = map(expr_Λ, values)
        vecs = reshape(vectors, χr, χl, Nk)
        v0 = reshape(init_vector, χr, χl)

        Onel = ones(χl, χl); Onel = convert(S, Onel)
        Oner = ones(χr, χr); Oner = convert(S, Oner)
        ∂y_∂Ql += -t * ein"n,n,kbn,ka,kk -> ab"(Λ, weight, conj(vecs), v0, Oner)
        ∂y_∂Qr += -t * ein"n,n,bkn,ak,kk -> ab"(Λ, weight, conj(vecs), v0, Onel)

        if U <: AbstractMatrix
            ∂y_∂Rl += -t * ein"n,n,kl,lbn,ka,kl -> ab"(Λ, weight, ψr.R, conj(vecs), v0, Oner)
            ∂y_∂Rr += -t * ein"n,n,kl,bln,ak,kl -> ab"(Λ, weight, ψl.R, conj(vecs), v0, Onel)
        else
            Onel = ones(χl, χl, size(ψ.R,3)); Onel = convert(U, Onel)
            Oner = ones(χr, χr, size(ψ.R,3)); Oner = convert(U, Oner)
            ∂y_∂Rl += -t * ein"n,n,klm,lbn,ka,klm -> ab"(Λ, weight, ψr.R, conj(vecs), v0, Oner)
            ∂y_∂Rr += -t * ein"n,n,klm,bln,ak,klm -> ab"(Λ, weight, ψl.R, conj(vecs), v0, Onel)
        end
    end
    
    factor = Ns/Nr
    y, ∂y_∂t = map(x -> x * factor, res)
    ∂y_∂Ql, ∂y_∂Qr, ∂y_∂Rl, ∂y_∂Rr = map(x -> x * factor, (∂y_∂Ql, ∂y_∂Qr, ∂y_∂Rl, ∂y_∂Rr))

    function logtrexp_pullback(ȳ)
        ∂y_∂t, ∂y_∂Ql, ∂y_∂Qr, ∂y_∂Rl, ∂y_∂Rr = map(x -> x/y, (∂y_∂t, ∂y_∂Ql, ∂y_∂Qr, ∂y_∂Rl, ∂y_∂Rr))

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