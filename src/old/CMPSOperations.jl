using Test, JuliaCMPO
using LinearAlgebra, FiniteTLanczos
using Random; Random.seed!()


estimator1 = nothing
estimator2 = TraceEstimator(full_ed_ftrace, FTLMOptions(FTLMOptions(), processor = processor))

options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=100, processor = processor)
estimator3 = TraceEstimator(standard_ftlm_ftrace, options)

estimator_list = [estimator1, estimator2]
β = rand()

@testset "normalize a CMPS" begin
    for vD = 1:2
        ψ = init_cmps(10, vD)
        for trace_estimator in estimator_list
            ψ1 = solver(x-> normalize(x, β, trace_estimator), ψ)
            @test solver(x->norm(x,β, trace_estimator), ψ1) ≈ 1.
        end

        ψ1 = solver(x-> normalize(x, β, estimator3), ψ)
        try
            solver(x->norm(x,β, estimator3), ψ1) 
            @test true
        catch
            @test false
        end
    end
end

@testset "logfidelity" begin
    for vD = 1:2
        ψ0 = init_cmps(10, vD)
        ψ1 = init_cmps(12, vD)
        fidel0 = solver((x1,x0)->logfidelity(x1, x0, β, nothing), ψ1, ψ0)
        for trace_estimator in estimator_list
            fidel1 = solver((x1,x0)->logfidelity(x1, x0, β, trace_estimator), ψ1, ψ0)
            @test ≈(fidel1, fidel0)
        end

        try
            solver((x1,x0)->logfidelity(x1, x0, β, estimator3), ψ1, ψ0)
            @test true
        catch
            @test false
        end
    end
end




