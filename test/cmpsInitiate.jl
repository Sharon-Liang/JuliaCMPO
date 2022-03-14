using Test, cMPO 

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


@testset "cmps_compress" begin
    β = 1.0
    χ = rand(6:10)
    vir_dim = rand(1:2)

    ψ0 = init_cmps(χ, vir_dim)
    ψ = cmps_compress(ψ0, χ-2, β)

    @test size(ψ.Q) == (χ-2, χ-2)
    @test size(ψ.R) == (χ-2, χ-2, vir_dim) 
    @test fidelity(ψ, ψ0, β) ≈ fidelity(ψ0, ψ0, β)
end

@testset "Initiate via boundary cmps" begin
    m = TFIsing(1.0,1.0)
    β = 1
    χ = rand([2,4,8,16])
    ψ = init_cmps(χ, m, β)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ)

    χ = rand(10:20)
    ψ = init_cmps(χ, m, β)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ) 
end