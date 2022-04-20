using Test, cMPO 
using Random; Random.seed!()
#solver = cpu_solver

@testset "Initiate a CMPS randomly" begin
    χ = rand(2:12)
    vir_dim = rand(2:4)

    ψ = init_cmps(χ)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ)

    ψ = init_cmps(χ, vir_dim)
    @test size(ψ.Q) == (χ, χ)
    @test size(ψ.R) == (χ, χ, vir_dim)
end

for solver in [cpu_solver, gpu_solver] 
    @testset "$(solver)" begin
        @testset "Initiate via boundary CMPS" begin
            β = 2.0
            n = 3
            for m in [TFIsing(1.0,1.0), TFIsing_2D_helical(1.0,1.0,2)], χ in [2^n, 2^n-3]
                solver == cpu_solver ? type = CMPS : type = CuCMPS
                ψ = solver(Tm->init_cmps(χ, Tm, β), m.Tmatrix)
                @test typeof(ψ) <: type
                @test size(ψ.Q) == (χ, χ)
                @test size(ψ.R)[1:2] == (χ, χ) && length(size(ψ.R)) == length(size(m.Tmatrix.R))
            end
        end
    end
end



