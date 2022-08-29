using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()

@testset "matrix logtrexp" begin
    D = 4; T = ComplexF64
    t = rand()
    for A in [randn(D,D), randn(T,D,D)]
        A = symmetrize(A)
        log_trexp = exp(t*A) |> tr |> real |> log
        @test logtrexp(t, solver(A)) ≈ log_trexp
    end
end

@testset "CMPSMatrix logtrexp" begin
    Nl, Nr = 10, 11
    ψl, ψr = init_cmps(Nl), init_cmps(Nr)
    Cmatrix = ψl * ψr |> solver
    M = Cmatrix |> tomatrix
    t = rand()


    log_trexp = logtrexp(t, M)
    @testset "Explicit construct CMPSMatrix" begin
        @test logtrexp(t, Cmatrix) ≈ log_trexp
    end

    @testset "Full ED method" begin
        trace_estimator = TraceEstimator(full_ed_ftrace, FTLMOptions(FTLMOptions(), processor = processor))
        res = logtrexp(t, Cmatrix, trace_estimator)
        @test res ≈ log_trexp
    end

    @testset "standard FTLM method" begin
        Ns = size(Cmatrix, 1)
        options = FTLMOptions(FTLMOptions(), processor = processor)
        trace_estimator = TraceEstimator(standard_ftlm_ftrace, options)
        res = logtrexp(t, Cmatrix, trace_estimator)
        @unpack Nr = options
        @test ≈(res, log_trexp, rtol=1/√Nr)
    end
    
    @testset "orthogonalized FTLM method" begin
        Ns = size(Cmatrix,1)
        options = FTLMOptions(FTLMOptions(), processor = processor)
        trace_estimator = TraceEstimator(orthogonalized_ftlm_ftrace, options)
        res = logtrexp(t, Cmatrix, trace_estimator)
        @unpack Nr = options
        @test ≈(res, log_trexp, rtol=1/Nr)
    end
end