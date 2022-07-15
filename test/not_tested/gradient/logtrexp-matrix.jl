#logtrexp
using Test, JuliaCMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()
using FiniteDifferences, Zygote

ngradient(f, x) = grad(central_fdm(5, 1), f, x)
zgradient = Zygote.gradient

processor = CPU
solver = solver_function(processor)

Nl, Nr = 3, 4

ψl, ψr = init_cmps(Nl), init_cmps(Nr)
Cmatrix = ψl * ψr
M = Cmatrix |> Matrix
t = rand()

M = solver(x->x, M)
Cmatrix = solver(x->x, Cmatrix)

∂t_M = x->logtrexp(x, M)
∂t_Cmatrix = x->logtrexp(x, Cmatrix)
∂Matrix = x->logtrexp(t, x)

∂t_Cmatrix_estimate = (x, trace_estimator)->logtrexp(x, Cmatrix, trace_estimator)
∂Matrix_estimate = (x, trace_estimator)->logtrexp(t, x, trace_estimator)

#Benchmarks
∂t_M_ndiff = ngradient(∂t_M, t)[1] 
∂Matrix_Cmatrix_ndiff = ngradient(∂Matrix, Cmatrix)[1]

@testset "logtrexp-matrix" begin
    ∂t_M_zygote = zgradient(∂t_M, t)[1]
    @test ≈(∂t_M_zygote, ∂t_M_ndiff; rtol = 1e-5, atol = 1e-5)

    ∂Matrix_M_ndiff = ngradient(∂Matrix, M)[1] 
    ∂Matrix_M_zygote = zgradient(∂Matrix, M)[1]
    @test ≈(∂Matrix_M_zygote, ∂Matrix_M_ndiff; rtol = 1e-5, atol = 1e-5)
end

@testset "logtrexp-CMPSMatrix" begin 
    ∂t_Cmatrix_ndiff = ngradient(∂t_Cmatrix, t)[1]  
    ∂t_Cmatrix_zygote = zgradient(∂t_Cmatrix, t)[1] 
    @test ≈(∂t_Cmatrix_ndiff, ∂t_M_ndiff; rtol = 1e-5, atol = 1e-5)
    @test ≈(∂t_Cmatrix_zygote, ∂t_Cmatrix_ndiff; rtol = 1e-5, atol = 1e-5)

    ∂Matrix_Cmatrix_zygote = zgradient(∂Matrix, Cmatrix)[1]
    @test ≈(∂Matrix_Cmatrix_zygote, ∂Matrix_Cmatrix_ndiff; rtol = 1e-5, atol = 1e-5)
end


@testset "logtrexp-CMPSMatrix-nothing" begin
    trace_estimator = nothing
    ∂t_Cmatrix_estimate_ndiff = ngradient(x->∂t_Cmatrix_estimate(x,trace_estimator), t)[1]  
    ∂t_Cmatrix_estimate_zygote = zgradient(x->∂t_Cmatrix_estimate(x,trace_estimator), t)[1] 
    @test ≈(∂t_Cmatrix_estimate_ndiff, ∂t_M_ndiff; rtol = 1e-5, atol = 1e-5)
    @test ≈(∂t_Cmatrix_estimate_zygote, ∂t_Cmatrix_estimate_ndiff; rtol = 1e-5, atol = 1e-5)


    ∂Matrix_Cmatrix_estimate_ndiff = ngradient(x->∂Matrix_estimate(x, trace_estimator), Cmatrix)[1]
    ∂Matrix_Cmatrix_estimate_zygote = zgradient(x->∂Matrix_estimate(x, trace_estimator), Cmatrix)[1]
    @test ≈(∂Matrix_Cmatrix_estimate_ndiff, ∂Matrix_Cmatrix_ndiff; rtol = 1e-5, atol = 1e-5)
    @test ≈(∂Matrix_Cmatrix_estimate_zygote, ∂Matrix_Cmatrix_estimate_ndiff; rtol = 1e-5, atol = 1e-5)
end


@testset "logtrexp-CMPSMatrix-Full_ED" begin
    trace_estimator = TraceEstimator(Full_ED, FTLMOptions(FTLMOptions(), processor = processor))
    ∂t_Cmatrix_estimate_zygote = zgradient(x->∂t_Cmatrix_estimate(x,trace_estimator), t)[1] 
    @test ≈(∂t_Cmatrix_estimate_zygote, ∂t_M_ndiff; rtol = 1e-5, atol = 1e-5)

    ∂Matrix_Cmatrix_estimate_zygote = zgradient(x->∂Matrix_estimate(x, trace_estimator), Cmatrix)[1]
    @test ≈(∂Matrix_Cmatrix_estimate_zygote, ∂Matrix_Cmatrix_ndiff; rtol = 1e-5, atol = 1e-5)
end


@testset "logtrexp-CMPSMatrix-FullSampling_FTLM" begin
    Ns = size(Cmatrix,1)
    options = FTLMOptions(FTLMOptions(), Ne = 0, Nk=Ns, processor = processor)
    trace_estimator = TraceEstimator(FullSampling_FTLM, options)
    ∂t_Cmatrix_estimate_zygote = zgradient(x->∂t_Cmatrix_estimate(x,trace_estimator), t)[1] 
    @test ≈(∂t_Cmatrix_estimate_zygote, ∂t_M_ndiff; rtol = 1e-5, atol = 1e-5)
    
    ∂Matrix_Cmatrix_estimate_zygote = zgradient(x->∂Matrix_estimate(x, trace_estimator), Cmatrix)[1]
    @test ≈(∂Matrix_Cmatrix_estimate_zygote, ∂Matrix_Cmatrix_ndiff; rtol = 1e-5, atol = 1e-5)
end

@testset "logtrexp-CMPSMatrix-simple_FTLM" begin
    Ns = size(Cmatrix,1)
    options = FTLMOptions(FTLMOptions(), processor = processor)
    trace_estimator = TraceEstimator(simple_FTLM, options)
    try
        ∂t_Cmatrix_estimate_zygote = zgradient(x->∂t_Cmatrix_estimate(x,trace_estimator), t)[1] 
        ∂Matrix_Cmatrix_estimate_zygote = zgradient(x->∂Matrix_estimate(x, trace_estimator), Cmatrix)[1]
        @test true
    catch
        @test false
    end
end


@testset "logtrexp-CMPSMatrix-orthogonalized_FTLM" begin
    Ns = size(Cmatrix,1)
    options = FTLMOptions(FTLMOptions(), Ne = Ns , processor = processor)
    trace_estimator = TraceEstimator(orthogonalized_FTLM, options)
    ∂t_Cmatrix_estimate_zygote = zgradient(x->∂t_Cmatrix_estimate(x,trace_estimator), t)[1] 
    @test ≈(∂t_Cmatrix_estimate_zygote, ∂t_M_ndiff; rtol = 1e-5, atol = 1e-5)
    
    ∂Matrix_Cmatrix_estimate_zygote = zgradient(x->∂Matrix_estimate(x, trace_estimator), Cmatrix)[1]
    @test ≈(∂Matrix_Cmatrix_estimate_zygote, ∂Matrix_Cmatrix_ndiff; rtol = 1e-5, atol = 1e-5)

    options = FTLMOptions(FTLMOptions(), processor = processor)
    trace_estimator = TraceEstimator(orthogonalized_FTLM, options)
    try
        ∂t_Cmatrix_estimate_zygote = zgradient(x->∂t_Cmatrix_estimate(x,trace_estimator), t)[1] 
        ∂Matrix_Cmatrix_estimate_zygote = zgradient(x->∂Matrix_estimate(x, trace_estimator), Cmatrix)[1]
        @test true
    catch
        @test false
    end
end

@testset "logtrexp-CMPSMatrix-replaced_FTLM" begin
    Ns = size(Cmatrix,1)
    options = FTLMOptions(FTLMOptions(), Ne = Ns, processor = processor)
    trace_estimator = TraceEstimator(replaced_FTLM, options)
    try
        ∂t_Cmatrix_estimate_zygote = zgradient(x->∂t_Cmatrix_estimate(x,trace_estimator), t)[1] 
        ∂Matrix_Cmatrix_estimate_zygote = zgradient(x->∂Matrix_estimate(x, trace_estimator), Cmatrix)[1]
        @test true
    catch
        @test false
    end
end