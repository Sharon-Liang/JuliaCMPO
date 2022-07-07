using Test, JuliaCMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()

Nl, Nr = 10, 15

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

    for solver in [cpu_solver, gpu_solver] 
        @testset "$(solver)" begin
            @testset "Full Rank itFOLM generate eigenvalues and eigenvectors" begin
                evals, evecs = eigen(M)
                N = size(Cmatrix, 1)
                for distr in [Gaussian, Rademacher, rand]
                    v0 = random_unit_vector(N, distr)
                    v0 = solver(x->x, v0)
                    Cmatrix = solver(x->x, Cmatrix)
                    res = itFOLM(Cmatrix, init_vector = v0, ncv = N) |> eigensolver
                    @unpack values, vectors = res
                    @test Vector(values) ≈ evals
                    @test Matrix(vectors * diagm(values) * vectors') ≈ M
                end
            end
        end
    end
end


#=
@testset "logtrexp function" begin
    t = rand()
    res0 = logtrexp(t,M)
    for method in [Full_ED, FullSampling_FTLM]
        @test logtrexp(t, Cmatrix, method) ≈ res0
    end
end
=#

#TODO: test itFOLMe