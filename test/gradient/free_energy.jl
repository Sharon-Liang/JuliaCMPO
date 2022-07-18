#free energy
using Test, JuliaCMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()
using FiniteDifferences, Zygote

import JuliaCMPO:free_energy

#ngradient(f, x) = grad(central_fdm(5, 1), f, x)
#zgradient = Zygote.gradient

#processor = CPU
#solver = solver_function(processor)

estimator1 = nothing
estimator2 = TraceEstimator(Full_ED, FTLMOptions(FTLMOptions(), processor = processor))

options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=100, processor = processor)
estimator3 = TraceEstimator(FullSampling_FTLM, options)

estimator_list = [estimator1, estimator2, estimator3]

@testset "TFIsing model" begin
    β = 2.0; χ = 4
    m = TFIsing(1.0, 1.0)
    Tm = solver(x->x, m.Tmatrix)
    for trace_estimator in estimator_list
        ψ = init_cmps(χ) |> diagQ
        dQ = convert(Vector, diag(ψ.Q))
        R = convert(Array, ψ.R)
        pars = Zygote.Params([dQ, R])
        loss0() = free_energy(solver(CMPS_generate, diagm(dQ), R), Tm, β, nothing)
        loss1() = free_energy(solver(CMPS_generate, diagm(dQ), R), Tm, β, trace_estimator)
        p0, f0, g0! = optim_functions(loss0, pars)  
        p1, f1, g1! = optim_functions(loss1, pars)
        grad_ndiff = ngradient(f0, p0)[1]
        grad_zygote = similar(p1); g1!(grad_zygote, p1)
        @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
    end
end