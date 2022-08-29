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
estimator2 = TraceEstimator(full_ed_ftrace, FTLMOptions(FTLMOptions(), processor = processor))

options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=100, processor = processor)
estimator3 = TraceEstimator(standard_ftlm_ftrace, options)

@testset "TFIsing model" begin
    β = 2.0; χ = 4
    m = TFIsing(1.0, 1.0) |> solver
    
    ψ = init_cmps(χ) |> diagQ
    dQ = convert(Vector, diag(ψ.Q))
    R = convert(Array, ψ.R)
    pars = Zygote.Params([dQ, R])
    loss1() = free_energy(solver(cmps_generate(diagm(dQ), R)), m, β, estimator1)
    loss2() = free_energy(solver(cmps_generate(diagm(dQ), R)), m, β, estimator2)
    loss3() = free_energy(solver(cmps_generate(diagm(dQ), R)), m, β, estimator3)
    p1, f1, g1! = optim_functions(loss1, pars)
    p2, f2, g2! = optim_functions(loss2, pars)
    p3, f3, g3! = optim_functions(loss3, pars)

    grad_ndiff = ngradient(f1, p1)[1]
    grad_zygote1 = similar(p1); g1!(grad_zygote1, p1)
    @test ≈(grad_zygote1, grad_ndiff; rtol = 1e-5, atol = 1e-5)
        
    grad_zygote2 = similar(p2); g2!(grad_zygote2, p2)
    @test ≈(grad_zygote2, grad_ndiff; rtol = 1e-5, atol = 1e-5)

    try 
        g2!(similar(p3), p3)
        @test true
    catch
        @test false
    end
end
