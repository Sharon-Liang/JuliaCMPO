using cMPO, Test
using LinearAlgebra
using Random; Random.seed!()

@testset "utilities for hermitian matrices" begin
    D = 4; T = ComplexF64
    for A in [randn(D,D), randn(T,D,D)]
        log_trexp = symmetrize(A) |> exp |> tr |> real |> log
        for solver in [cpu_solver, gpu_solver]
            @test ishermitian(solver(symmetrize, A))
            @test solver(logtrexp, A) ≈ log_trexp
        end
    end
end

