using Test, JuliaCMPO 
using Random; Random.seed!()

@testset "Initiate a CMPS randomly" begin
    χ = rand(2:12)
    vir_dim = rand(2:4)

    ψ = init_cmps(χ)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ)

    ψ = init_cmps(χ, vir_dim)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ, vir_dim)
end

@testset "Initiate via boundary CMPS" begin
    estimator1 = nothing
    estimator2 = TraceEstimator(full_ed_ftrace, FTLMOptions(FTLMOptions(), processor = processor))

    estimator_list = [estimator1, estimator2]

    β = rand()
    n = 3
    for m in [TFIsing(1.0,1.0), TFIsing_2D_helical(1.0,1.0,2)], χ in [2^n, 2^n-3]
        processor == CPU ? tensor_type = CMPS : tensor_type = CuCMPS
        for trace_estimator in estimator_list
            compress_options = CompressOptions(;trace_estimator, processor)
            ψ = init_cmps(χ, solver(m), β; compress_options)
            @test typeof(ψ) <: tensor_type
            @test size(ψ.Q) == (χ, χ)
            @test size(ψ.R)[1:2] == (χ, χ) && length(size(ψ.R)) == length(size(m.R))
        end
    end
end



