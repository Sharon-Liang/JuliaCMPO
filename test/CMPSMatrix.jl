using Test, JuliaCMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()

Nlist = [10, 20]
Nl, Nr = rand(Nlist), rand(Nlist)

ψl, ψr = init_cmps(Nl), init_cmps(Nr)
Cmatrix = ψl * ψr
P0 = Cmatrix |> Matrix

@testset "size" begin
    @test size(P0) == size(Cmatrix)
end

@testset "getindex" begin
    @test all(map((i,j)->P0[i,j] ≈ getindex(Cmatrix,i,j), 1:Nl, 1:Nr))
    @test all(map(i->P0[i] ≈ getindex(Cmatrix,i), 1:length(P0)))
end

@testset "Full Rank itFOLM generate eigenvalues and eigenvectors" begin
    evals, evecs = eigen(P0)
    N = size(Cmatrix)[1]
    for dist in [Gaussian, Rademacher, rand]
        v0 = random_unit_vector(N, dist)
        res = itFOLM(Cmatrix, init_vector = v0, ncv = N) |> eigensolver
        @unpack values, vectors = res
        @test values ≈ evals
        @test vectors * diagm(values) * vectors' ≈ P0
    end
end

@testset "logtrexp function" begin
    t = rand()
    res0 = logtrexp(t,P0)
    for method in [Full_ED, FullSampling_FTLM]
        @test logtrexp(t, Cmatrix, method) ≈ res0
    end
end

#TODO: test itFOLMe