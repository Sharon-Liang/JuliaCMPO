using Test, JuliaCMPO
using LinearAlgebra, FiniteTLanczos
using Random; Random.seed!()


estimator1 = nothing
estimator2 = TraceEstimator(Full_ED, FTLMOptions(FTLMOptions(), processor = processor))

options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=100, processor = processor)
estimator3 = TraceEstimator(FullSampling_FTLM, options)

estimator_list = [estimator1, estimator2, estimator3]
β = rand()

@testset "normalize a CMPS" begin
    for vD = 1:2, trace_estimator in estimator_list
        ψ = init_cmps(10, vD)
        ψ1 = solver(x-> normalize(x, β, trace_estimator), ψ)
        @test solver(x->norm(x,β, trace_estimator), ψ1) ≈ 1.
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
    end
end




