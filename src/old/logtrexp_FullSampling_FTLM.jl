"""
    rrule for logtrexp(tM) function where `typeof(M)=CMPSMatrix`, 
    `trace_estimator = FullSampling_FTLM`
"""
function ChainRules.rrule(::typeof(logtrexp), 
                          t::Real, M::CMPSMatrix{Ts,T,S,U}, 
                          trace_estimator::TraceEstimator{Tf, To}
                          ) where {Ts,T,S,U,Tf<:typeof(FullSampling_FTLM),To}
    @unpack options = trace_estimator
    @unpack Nk, processor = options
    @unpack ψl, ψr = M
    Ns = size(M, 1)
    Nk = min(Nk, Ns)

    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)
    solver = solver_function(processor)
    
    sign(t) == 1 ? which = :LR : which=:SR
    processor == CPU ? x0 = rand(T,Ns) : x0 = CUDA.rand(T, Ns)
    e0, _, _ =eigsolve(M, x0, 1, which, ishermitian = true, tol=1.e-12)
    e1 = e0[1]
    expr_Λ = e -> exp(t*(e-e1))
    expr_∂y_∂t = e -> e * exp(t*(e-e1))
    expr = (expr_Λ, expr_∂y_∂t)

    ortho_basis = basis_generate(Ns)
    res = zeros(2)
    if processor == CPU
        ∂y_∂Ql, ∂y_∂Qr = zeros(T, size(ψl.Q)), zeros(T, size(ψr.Q))
        ∂y_∂Rl, ∂y_∂Rr = zeros(T, size(ψl.R)), zeros(T, size(ψr.R))
    else
        ∂y_∂Ql, ∂y_∂Qr = CUDA.zeros(T, size(ψl.Q)), CUDA.zeros(T, size(ψr.Q))
        ∂y_∂Rl, ∂y_∂Rr = CUDA.zeros(T, size(ψl.R)), CUDA.zeros(T, size(ψr.R))
    end

    for r = 1: Ns
        v0 = ortho_basis[:, r]
        v0 = solver(x->x, v0)
        @unpack init_vector, weight, values, vectors = itFOLM(M, init_vector = v0, Nk = Nk) |> eigensolver
        #Nk = size(values,1)
        
        func = f -> begin
            Λ = map(f, values)
            Z = ein"i,i,i -> "(Λ, weight, conj(weight))
            return Array(Z)[1]
        end
        res = map(+, res, map(func, expr))

        Λ = map(expr_Λ, values)
        vectors = reshape(vectors, χr, χl, Nk)
        v0 = reshape(init_vector, χr, χl)

        if processor == CPU
            Onel = ones(T, χl, χl)
            Oner = ones(T, χr, χr)
        else
            Onel = CUDA.ones(T, χl, χl)
            Oner = CUDA.ones(T, χr, χr)
        end
        ∂y_∂Ql_temp = -t * ein"n,n,kbn,ka,kk -> ab"(Λ, weight, conj(vectors), v0, Oner)
        ∂y_∂Qr_temp = -t * ein"n,n,bkn,ak,kk -> ab"(Λ, weight, conj(vectors), v0, Onel)
        ∂y_∂Ql = map(+, ∂y_∂Ql, ∂y_∂Ql_temp)
        ∂y_∂Qr = map(+, ∂y_∂Qr, ∂y_∂Qr_temp)

        if U <: AbstractMatrix
            ∂y_∂Rl_temp = -t * ein"n,n,kl,lbn,ka,kl -> ab"(Λ, weight, ψr.R, conj(vectors), v0, Oner)
            ∂y_∂Rr_temp = -t * ein"n,n,kl,bln,ak,kl -> ab"(Λ, weight, ψl.R, conj(vectors), v0, Onel)
        else
            if processor == CPU
                Onel = ones(T, χl, χl, size(ψl.R,3))
                Oner = ones(T, χr, χr, size(ψr.R,3))
            else
                Onel = CUDA.ones(T, χl, χl, size(ψl.R,3))
                Oner = CUDA.ones(T, χr, χr, size(ψr.R,3))
            end
            ∂y_∂Rl_temp = -t * ein"n,n,klm,lbn,ka,klm -> abm"(Λ, weight, ψr.R, conj(vectors), v0, Oner)
            ∂y_∂Rr_temp =  -t * ein"n,n,klm,bln,ak,klm -> abm"(Λ, weight, ψl.R, conj(vectors), v0, Onel)
        end
        ∂y_∂Rl = map(+, ∂y_∂Rl, ∂y_∂Rl_temp)
        ∂y_∂Rr = map(+, ∂y_∂Rr, ∂y_∂Rr_temp)
    end
    
    y, ∂y_∂t = res
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