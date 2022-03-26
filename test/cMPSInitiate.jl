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
    m = TFIsing(1.0,1.0)
    β = 2.0
    χ = rand(3:6)
    ψ = init_cmps(χ, m.Tmatrix, β)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ) 
end




