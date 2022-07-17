using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()

@testset "matrix logtrexp" begin
    D = 4; T = ComplexF64
    t = rand()
    for A in [randn(D,D), randn(T,D,D)]
        A = symmetrize(A)
        log_trexp = exp(t*A) |> tr |> real |> log
        @test solver(x->logtrexp(t,x), A) ≈ log_trexp
    end
end

@testset "CMPSMatrix logtrexp" begin
    Nl, Nr = 10, 11
    ψl, ψr = init_cmps(Nl), init_cmps(Nr)
    Cmatrix = ψl * ψr; Cmatrix = solver(x->x, Cmatrix)
    M = Cmatrix |> Matrix; M = solver(x->x, M)
    t = rand()


    log_trexp = logtrexp(t, M)
    @testset "Explicit construct CMPSMatrix" begin
        @test logtrexp(t, Cmatrix) ≈ log_trexp
    end

    @testset "Full_ED" begin
        trace_estimator = TraceEstimator(Full_ED, FTLMOptions(FTLMOptions(), processor = processor))
        res = logtrexp(t, Cmatrix, trace_estimator)
        @test res ≈ log_trexp
    end

    @testset "simple_FTLM method" begin
        Ns = size(Cmatrix,1)
        options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=Ns, processor = processor)
        trace_estimator = TraceEstimator(FullSampling_FTLM, options)
        res = logtrexp(t, Cmatrix, trace_estimator)
        @test res ≈ log_trexp

        options = FTLMOptions(FTLMOptions(), processor = processor)
        trace_estimator = TraceEstimator(simple_FTLM, options)
        res = logtrexp(t, Cmatrix, trace_estimator)
        @unpack Nr = options
        @test ≈(res, log_trexp, rtol=1/√Nr)
    end
    
    @testset "orthogonalized_FTLM method" begin
        Ns = size(Cmatrix,1)
        options = FTLMOptions(FTLMOptions(), processor = processor)
        trace_estimator = TraceEstimator(orthogonalized_FTLM, options)
        res = logtrexp(t, Cmatrix, trace_estimator)
        @unpack Nr = options
        @test ≈(res, log_trexp, rtol=1/Nr)
    end
end