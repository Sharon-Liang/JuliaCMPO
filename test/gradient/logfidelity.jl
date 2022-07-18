#logfidelity
using Test, JuliaCMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()
using FiniteDifferences, Zygote

#ngradient(f, x) = grad(central_fdm(5, 1), f, x)
#zgradient = Zygote.gradient

#processor = CPU
#solver = solver_function(processor)

estimator1 = nothing
estimator2 = TraceEstimator(Full_ED, FTLMOptions(FTLMOptions(), processor = processor))

options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=100, processor = processor)
estimator3 = TraceEstimator(FullSampling_FTLM, options)

estimator_list = [estimator1, estimator2, estimator3]

β = rand()
χ = 4
for vD = 1:2, trace_estimator in estimator_list
    ψ0 = solver(x->x, init_cmps(χ+2, vD))
    ψ1 = init_cmps(χ, vD) |> diagQ
    dQ = convert(Vector, diag(ψ1.Q))
    R = convert(Array, ψ1.R)
    pars = Zygote.Params([dQ, R])
    loss0() = -logfidelity(solver(CMPS_generate, diagm(dQ), R), ψ0, β, nothing)
    loss1() = -logfidelity(solver(CMPS_generate, diagm(dQ), R), ψ0, β, trace_estimator)
    p0, f0, g0! = optim_functions(loss0, pars)
    p1, f1, g1! = optim_functions(loss1, pars)

    grad_ndiff = ngradient(f0, p0)[1]
    grad_zygote = similar(p1); g1!(grad_zygote, p1)
    @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
end

