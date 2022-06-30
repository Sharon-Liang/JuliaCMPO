using Test, cMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()

Nlist = [10, 20]
Nl, Nr = rand(Nlist), rand(Nlist)

ψl, ψr = init_cmps(Nl), init_cmps(Nr)
Cmatrix = CMPSMatrix(ψl, ψr)

P0 = ψl * ψr

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
    for dist in [FiniteTLanczos.Gaussian, FiniteTLanczos.Rademacher, rand]
        v0 = random_unit_vector(N, dist)
        res = itFOLM(Cmatrix, init_vector = v0, ncv = N) |> eigensolver
        @unpack values, vectors = res
        @test values ≈ evals
        @test vectors * diagm(values) * vectors' ≈ P0
    end
end

@testset "logtrexp function" begin
    res0 = logtrexp(P0)
    for method in [FiniteTLanczos.Full_ED, FiniteTLanczos.FullSampling_FTLM]
        @test logtrexp(Cmatrix, method) ≈ res0
    end
end

#TODO: test itFOLMe