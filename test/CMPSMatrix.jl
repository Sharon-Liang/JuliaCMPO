using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()

Nl, Nr = 10, 11

for vD = 1:2
    ψl, ψr = init_cmps(Nl, vD), init_cmps(Nr, vD)
    Cmatrix = ψl * ψr 
    M = Cmatrix |> tomatrix 

    @testset "size" begin
        @test size(M) == size(Cmatrix)
    end

    @testset "getindex" begin
        @test all(map((i,j)->M[i,j] ≈ getindex(Cmatrix,i,j), 1:Nl, 1:Nr))
        @test all(map(i->M[i] ≈ getindex(Cmatrix,i), 1:length(M)))
    end

    @testset "*" begin
        v0 = rand(size(M,1)) |> solver
        m0 = rand(size(M,1), 3) |> solver

        @test ≈(solver(M) * v0, solver(Cmatrix) * v0)
        @test ≈(solver(M) * m0, solver(Cmatrix) * m0)
    end

    @testset "Full Rank Iterative Full Orthogonalized Lanczos Algorithm generate exact eigen pairs" begin
        evals, evecs = eigen(M)
        Nk = size(Cmatrix, 1)
        init_vector = random_unit_vector(Nk, rand) |> solver
        res = fullortho_lanczos(solver(Cmatrix); init_vector, Nk) |> eigensolver
        @unpack values, vectors = res
        @test Vector(values) ≈ evals
        @test Matrix(vectors * diagm(values) * vectors') ≈ M
    end
end


