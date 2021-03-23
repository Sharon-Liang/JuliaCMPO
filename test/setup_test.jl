using LinearAlgebra
using Random; Random.seed!()

@testset "utilities for hermitian matrices" begin
    A = rand(ComplexF32,(2,2))
    sa = symmetrize(A)
    tr_exp = exp(sa) |> tr |> real
    log_trexp = log(tr_exp)

    @test ishermitian(sa)
    @test isapprox(trexp(sa) |> value, tr_exp, rtol=1e-5)
    @test isapprox(logtrexp(sa), log_trexp, rtol=1e-5)
end

@testset "multiplications of cmps and cmpo" begin
    x = rand(2,2) |> symmetrize
    z = rand(2,2) |> symmetrize
    s = cmps(x, z)
    o = cmpo(x, x, z, zeros(2,2))

    s_arr = zeros(2,2,2)
    s_arr[:,:,1] = x; s_arr[:,:,2] = z

    i2 = Matrix(1.0I,2,2)
    ss = -(i2 ⊗ x + x ⊗ i2 + z ⊗ z)

    osq = i2 ⊗ x + x ⊗ i2 + z ⊗ z
    osr = x ⊗ i2

    soq = i2 ⊗ x + x ⊗ i2 + z ⊗ x
    sor = i2 ⊗ z

    ooq = i2 ⊗ x + x ⊗ i2 + z ⊗ x
    ool = i2 ⊗ z
    oor = x ⊗ i2
    oop = zeros(4,4)

    @test isapprox(s_arr, toarray(s))
    @test isapprox(s*s, ss)
    @test isapprox((o*s).Q, osq)
    @test isapprox((o*s).R, osr)
    @test isapprox((s*o).Q, soq)
    @test isapprox((s*o).R, sor)
    @test isapprox((o*o).Q, ooq)
    @test isapprox((o*o).R, oor)
    @test isapprox((o*o).L, ool)
    @test isapprox((o*o).P, oop)
end

@testset "normalize" begin
    s = init_cmps(2); β = 20
    ns = normalize(s, β)
    @test isapprox(ovlp(ns, β), 1, rtol=1.e-5)
end
