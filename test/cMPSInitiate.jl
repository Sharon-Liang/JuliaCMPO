using Test, cMPO 
using Random; Random.seed!()

@testset "Initiate a cMPS randomly" begin
    χ = rand(2:12)
    vir_dim = rand(2:4)

    ψ = init_cmps(χ)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ)

    ψ = init_cmps(χ, vir_dim)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ, vir_dim)
end


@testset "Initiate via boundary cmps" begin
    β = 2.0
    for m in [TFIsing(1.0,1.0), TFIsing_2D_helical(1.0,1.0,2)], χ in [4, 5],
        solver in [cpu_solver, gpu_solver]
        solver == cpu_solver ? type = CMPS : type = CuCMPS
        ψ = solver(x->init_cmps(χ, x, β), m.Tmatrix)
        @test typeof(ψ) <: type
        @test size(ψ.Q) == (χ, χ)
        @test size(ψ.R)[1:2] == (χ, χ) && length(size(ψ.R)) == length(size(m.Tmatrix.R))
    end
end



