using cMPO, Test
using LinearAlgebra
using Random; Random.seed!()

@testset "utilities for hermitian matrices" begin
    D = 4; T = ComplexF64
    for A in [randn(D,D), randn(T,D,D)], device in [:cpu, :gpu]
        sa = symmetrize(A)
        tr_exp = exp(sa) |> tr |> real
        log_trexp = log(tr_exp)

        @test ishermitian(sa)
        @test logtrexp(A, device=device) ≈ log_trexp
    end
end


@testset "-β * eigvals(A) |> sum and eigvals(-β * A) |> sum" begin
    D = 4; T = ComplexF64
    β = 10
    for A in [randn(D,D), randn(T,D,D)]
        A = A |> symmetrize |> Hermitian
        e1 = eigvals(A)
        e2 = eigvals(-β * A)
        a1 = exp.(-β .* e1) |> sum
        a2 = exp.(e2) |> sum
        @test a1 ≈ a2
    end
end

