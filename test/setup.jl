using cMPO, Test
using LinearAlgebra
using Random; Random.seed!()

@testset "utilities for hermitian matrices" begin
    D = 4; T = ComplexF64
    for A in [randn(D,D), randn(T,D,D)]
        sa = symmetrize(A)
        tr_exp = exp(sa) |> tr |> real
        log_trexp = log(tr_exp)

        @test ishermitian(sa)
        @test isapprox(trexp(sa) |> value, tr_exp, rtol=1e-5)
        @test isapprox(logtrexp(sa), log_trexp, rtol=1e-5)
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

@testset "normalize" begin
    T = ComplexF64
    D = 4
    for T in [Float64, ComplexF64]
        s = init_cmps(D, dtype=T); β = 20
        ns = normalize(s, β)
        @test isapprox(ovlp(ns, β), 1, rtol=1.e-5)
    end
end
