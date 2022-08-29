#project function
using Test, JuliaCMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()
using FiniteDifferences, Zygote

#ngradient(f, x) = grad(central_fdm(5, 1), f, x)
#zgradient = Zygote.gradient

#processor = CPU
#solver = solver_function(processor)

estimator1 = nothing
estimator2 = TraceEstimator(full_ed_ftrace, FTLMOptions(FTLMOptions(), processor = processor))

options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=100, processor = processor)
estimator3 = TraceEstimator(standard_ftlm_ftrace, options)

estimator_list = [estimator1, estimator2]


β = rand()
χ = 4
D = 2
loss0(v, ψ) = logfidelity(project(ψ, v), ψ, β, nothing)
          
@testset "project function" begin 
    for vD = 1:2
        processor == CPU ? u = rand(χ, D) : u = CUDA.rand(Float64, χ, D)
        ψ0 = init_cmps(χ, vD) |> solver
        grad_ndiff = ngradient(u->loss0(Array(u), CTensor(ψ0)), Array(u))[1]  
        for trace_estimator in estimator_list
            loss(v, ψ) = logfidelity(project(ψ, v), ψ, β, trace_estimator)
            grad_zygote = zgradient(u->loss(u, ψ0), u)[1] |> Array
            @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
        end
        try
            loss(v, ψ) = logfidelity(project(ψ, v), ψ, β, estimator3)
            grad_zygote = zgradient(u->loss(u, ψ0), u)[1] |> Array
            @test true
        catch
            @test false
        end

    end
end