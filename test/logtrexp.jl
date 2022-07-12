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
    Nl, Nr = 16, 17
    ψl, ψr = init_cmps(Nl), init_cmps(Nr)
    Cmatrix = ψl * ψr; Cmatrix = solver(x->x, Cmatrix)
    M = Cmatrix |> Matrix; M = solver(x->x, M)
    t = rand()


    log_trexp = logtrexp(t, M)
    @testset "Explicit construct CMPSMatrix" begin
        @test logtrexp(t, Cmatrix) ≈ log_trexp
    end

    @testset "Full_ED" begin
        estimator = TraceEstimator(Full_ED, FTLMOptions(FTLMOptions(), device = device))
        res = logtrexp(t, Cmatrix, estimator)
        @test res ≈ log_trexp
    end

    @testset "simple_FTLM method" begin
        Ns = size(Cmatrix,1)
        options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=Ns, device = device)
        estimator = TraceEstimator(FullSampling_FTLM, options)
        res = logtrexp(t, Cmatrix, estimator)
        @test res ≈ log_trexp

        options = FTLMOptions(FTLMOptions(), device = device)
        estimator = TraceEstimator(simple_FTLM, options)
        res = logtrexp(t, Cmatrix, estimator)
        @unpack Nr = options
        @test ≈(res, log_trexp, rtol=1/√Nr)
    end
    
    @testset "orthogonalized_FTLM method" begin
        Ns = size(Cmatrix,1)
        options = FTLMOptions(FTLMOptions(), device = device)
        estimator = TraceEstimator(orthogonalized_FTLM, options)
        res = logtrexp(t, Cmatrix, estimator)
        @unpack Nr = options
        @test ≈(res, log_trexp, rtol=1/Nr)
    end
end