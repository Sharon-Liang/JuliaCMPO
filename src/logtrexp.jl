"""
    logtrexp(t, M[, estimator::Union{Nothingm=, TraceEstimator}])

Calculate ``lnTr(e^{tM})``, where ``t`` is a real number and 
``M`` is matrix-like instances and is required to be hermitian.
"""
function logtrexp end
    
function logtrexp(t::Real, M::AbstractMatrix)
    M = symmetrize(M)
    vals, _ = eigensolver(M)
    return logsumexp(t*vals)
end
logtrexp(t::Real, M::CMPSMatrix) = logtrexp(t, tomatrix(M))
logtrexp(t, M, ::Nothing) = logtrexp(t, M)


function logtrexp(t::Real, M, trace_estimator::TraceEstimator)
    sign(t) == 1 ? which = :LR : which=:SR
    @unpack estimator, options = trace_estimator
    @unpack processor = options
    solver = solver_functions(processor)
    M = solver(M)

    processor == CPU ? x0 = rand(size(M,1)) : x0 = CUDA.rand(Float64, size(M,1))
    e0, _, _ =eigsolve(M, x0, 1, which, ishermitian = true, tol=1.e-12)
    e1 = e0[1] #typeof(e0) = Vector{Float64}
    expr = e -> exp(t * (e - e1))

    options = FTLMOptions(options, which = which)
    res = FiniteTLanczos.evaluate(TraceEstimator(estimator, options), M, expr)
    return log(res) + t*e1
end

