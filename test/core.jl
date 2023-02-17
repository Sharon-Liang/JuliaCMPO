using Test, JuliaCMPO 
using Random; Random.seed!()

#=
### *CMPS Initiate* 
=#
@testset "Initiate a CMPS randomly" begin
    χ = rand(2:12)
    D = rand(2:4)

    ψ = init_cmps(χ)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ)

    ψ = init_cmps(χ, D)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ, D)
end


@testset "Initiate via boundary CMPS" begin
    β = rand()
    n = 3
    for m in [model(TFIsingChain(), 1.0), model(TFIsingSquareHelical(), 1.0, 2)], χ in [2^n, 2^n-3]
    
        ψ = init_cmps(χ, solver(m), β)
        @test typeof(ψ.Q) <: typeof(m.Q)
        @test size(ψ.Q) == (χ, χ)
        @test size(ψ.R)[1:2] == (χ, χ) && length(size(ψ.R)) == length(size(m.R))
    end
end



