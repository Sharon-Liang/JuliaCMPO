"""
    rrule for logtrexp(tM) function where `typeof(M)=CMPSMatrix`, 
    `trace_estimator = orthogonalized_FTLM`
"""
function ChainRules.rrule(::typeof(logtrexp), 
                          t::Real, M::CMPSMatrix{Ts,T,S,U}, 
                          trace_estimator::TraceEstimator{Tf, To}
                          ) where {Ts,T,S,U,Tf<:typeof(orthogonalized_FTLM),To}
    @unpack options = trace_estimator
    @unpack distr, Nr, Nk, Ne, processor = options
    @unpack ψl, ψr = M
    solver = solver_function(processor)

    Ns = size(M, 1)
    Ne = min(Ne, Ns)
    χl, χr = size(ψl.Q, 1), size(ψr.Q, 1)

    sign(t) == 1 ? which = :LR : which=:SR
    processor == CPU ? x0 = rand(T,Ns) : x0 = CUDA.rand(T, Ns)
    e0, _, _ = eigsolve(M, x0, 1, which, ishermitian = true)
    e1 = e0[1]
    expr_Λ = e -> exp(t*(e-e1))
    expr_∂y_∂t = e -> e * exp(t*(e-e1))
    expr = (expr_Λ, expr_∂y_∂t)

    #exact eigen pair part
    krylovdim = max(30, Ne+1)
    vals, vecs, _ = eigsolve(M, x0, Ne, which, ishermitian = true, krylovdim=krylovdim)
    eigen_vals = solver(x->x, vals[1:Ne])
    eigen_vecs = hcat(vecs[1:Ne]...)

    func1 = f -> reduce(+, map(f, eigen_vals))
    res1 = map(func1, expr) #store y1, ∂y_∂t1
    
    Λ1 = map(expr_Λ, eigen_vals)
    eigen_vecs2 = reshape(eigen_vecs, χr, χl, Ne)
    Onel = ones(χl, χl); Onel = convert(S, Onel)
    Oner = ones(χr, χr); Oner = convert(S, Oner)

    ∂y_∂Ql1 = -t * ein"n,kbn,kan,kk -> ab"(Λ1, conj(eigen_vecs2), eigen_vecs2, Oner)
    ∂y_∂Qr1 = -t * ein"n,bkn,akn,kk -> ab"(Λ1, conj(eigen_vecs2), eigen_vecs2, Onel)
    if U <: AbstractMatrix
        ∂y_∂Rl1 = -t * ein"n,kl,lbn,kan,kl -> ab"(Λ1, ψr.R, conj(eigen_vecs2), eigen_vecs2, Oner)
        ∂y_∂Rr1 = -t * ein"n,kl,bln,akn,kl -> ab"(Λ1, ψl.R, conj(eigen_vecs2), eigen_vecs2, Onel)
    else
        Onel = ones(χl, χl, size(ψ.R,3)); Onel = convert(U, Onel)
        Oner = ones(χr, χr, size(ψ.R,3)); Oner = convert(U, Oner)
        ∂y_∂Rl1 = -t * ein"n,klm,lbn,kan,kl -> ab"(Λ1, ψr.R, conj(eigen_vecs2), eigen_vecs2, Oner)
        ∂y_∂Rr1 = -t * ein"n,klm,bln,akn,kl -> ab"(Λ1, ψl.R, conj(eigen_vecs2), eigen_vecs2, Onel)
    end


    #lanczos part
    res2 = zeros(2)
    ∂y_∂Ql2, ∂y_∂Qr2 = zeros(size(ψl.Q)), zeros(size(ψr.Q))
    ∂y_∂Rl2, ∂y_∂Rr2 = zeros(size(ψl.R)), zeros(size(ψr.R))
    for r = 1: Nr
        v0 = random_unit_vector(Ns, distr)
        v0 = solver(x->x, v0)
        @unpack init_vector, weight, values, vectors = 
            itFOLM(M, eigen_vecs, init_vector = v0, Nk = Nk) |> eigensolver
        Nk = size(values,1)
        
        func2 = f -> begin
            Λ = map(f, values)
            Z = ein"i,i,i -> "(Λ, weight, conj(weight))
            return Array(Z)[1]
        end
        res2 = map(+, res2, map(func2, expr))

        Λ2 = map(expr_Λ, values)
        vectors = reshape(vectors, χr, χl, Nk)
        v0 = reshape(init_vector, χr, χl)

        Onel = ones(χl, χl); Onel = convert(S, Onel)
        Oner = ones(χr, χr); Oner = convert(S, Oner)
        ∂y_∂Ql2_temp = -t * ein"n,n,kbn,ka,kk -> ab"(Λ2, weight, conj(vectors), v0, Oner)
        ∂y_∂Qr2_temp = -t * ein"n,n,bkn,ak,kk -> ab"(Λ2, weight, conj(vectors), v0, Onel)
        ∂y_∂Ql2 = map(+, ∂y_∂Ql2, ∂y_∂Ql2_temp)
        ∂y_∂Qr2 = map(+, ∂y_∂Qr2, ∂y_∂Qr2_temp)

        if U <: AbstractMatrix
            ∂y_∂Rl2_temp = -t * ein"n,n,kl,lbn,ka,kl -> ab"(Λ2, weight, ψr.R, conj(vectors), v0, Oner)
            ∂y_∂Rr2_temp = -t * ein"n,n,kl,bln,ak,kl -> ab"(Λ2, weight, ψl.R, conj(vectors), v0, Onel)
        else
            Onel = ones(χl, χl, size(ψ.R,3)); Onel = convert(U, Onel)
            Oner = ones(χr, χr, size(ψ.R,3)); Oner = convert(U, Oner)
            ∂y_∂Rl2_temp = -t * ein"n,n,klm,lbn,ka,klm -> ab"(Λ2, weight, ψr.R, conj(vectors), v0, Oner)
            ∂y_∂Rr2_temp = -t * ein"n,n,klm,bln,ak,klm -> ab"(Λ2, weight, ψl.R, conj(vectors), v0, Onel)
        end
        ∂y_∂Rl2 = map(+, ∂y_∂Rl2, ∂y_∂Rl2_temp)
        ∂y_∂Rr2 = map(+, ∂y_∂Rr2, ∂y_∂Rr2_temp)
    end

    factor = (Ns-Ne)/Nr
    res2 = map(x -> x * factor, res2)
    ∂y_∂Ql2, ∂y_∂Qr2, ∂y_∂Rl2, ∂y_∂Rr2 = map(x -> x * factor, (∂y_∂Ql2, ∂y_∂Qr2, ∂y_∂Rl2, ∂y_∂Rr2))

    y, ∂y_∂t = map(+, res1, res2)
    ∂y_∂Ql, ∂y_∂Qr, ∂y_∂Rl, ∂y_∂Rr = map(+,(∂y_∂Ql1, ∂y_∂Qr1, ∂y_∂Rl1, ∂y_∂Rr1),
                                            (∂y_∂Ql2, ∂y_∂Qr2, ∂y_∂Rl2, ∂y_∂Rr2))
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