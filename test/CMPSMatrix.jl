using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()

Nl, Nr = 10, 11

for vD = 1:2
    ψl, ψr = init_cmps(Nl, vD), init_cmps(Nr, vD)
    Cmatrix = ψl * ψr
    M = Cmatrix |> Matrix

    @testset "size" begin
        @test size(M) == size(Cmatrix)
    end

    @testset "getindex" begin
        @test all(map((i,j)->M[i,j] ≈ getindex(Cmatrix,i,j), 1:Nl, 1:Nr))
        @test all(map(i->M[i] ≈ getindex(Cmatrix,i), 1:length(M)))
    end

    @testset "CMPSMatrix: Full Rank itFOLM generate exact eigen pairs" begin
        evals, evecs = eigen(M)
        N = size(Cmatrix, 1)
        v0 = random_unit_vector(N); v0 = solver(x->x, v0)
        Cmatrix = solver(x->x, Cmatrix)
        res = itFOLM(Cmatrix, init_vector = v0, Nk = N) |> eigensolver
        @unpack values, vectors = res
        @test Vector(values) ≈ evals
        @test Matrix(vectors * diagm(values) * vectors') ≈ M
    end
end


