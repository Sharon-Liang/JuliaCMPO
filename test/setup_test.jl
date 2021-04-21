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

@testset "-β * eigvals(A) |> sum and eigvals(-β * A) |> sum" begin
    A = rand(4,4) |> symmetrize |> Hermitian
    β = 10
    e1 = eigvals(A)
    e2 = eigvals(-β * A)
    a1 = exp.(-β .* e1) |> sum
    a2 = exp.(e2) |> sum
    @test isapprox(a1, a2)
end

@testset "multiplications of cmps and cmpo: D-1 = 1" begin
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

    sos = -(i2 ⊗ osq + x ⊗ i2 ⊗ i2 + z ⊗ osr)

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
    @test isapprox(s*o*s, sos)
end


@testset "normalize" begin
    s = init_cmps(2); β = 20
    ns = normalize(s, β)
    @test isapprox(ovlp(ns, β), 1, rtol=1.e-5)
end

@testset "multiplications of cmps and cmpo: D-1 > 1" begin
    dtype = ComplexF64
    x = rand(dtype,(2,2)) |> symmetrize
    z = rand(2,2) |> symmetrize
    R = zeros(dtype,(2,2,2)); R[:,:,1] = x ; R[:,:,2] = x
    L = zeros(dtype,(2,2,2)); L[:,:,1] = z ; L[:,:,2] = z
    s = cmps(x, R)
    o = cmpo(x, R, L, zeros(dtype,(2,2,2,2)))

    s_arr = zeros(dtype,(2,2,3))
    s_arr[:,:,1] = x; s_arr[:,:,2] = x; s_arr[:,:,3] = x;

    i2 = Matrix(1.0I,2,2)
    ss = -(i2 ⊗ x + x ⊗ i2 + x ⊗ x + x ⊗ x)

    osq = i2 ⊗ x + x ⊗ i2 + z ⊗ x + z ⊗ x
    osr = zeros(dtype,(4,4,2))
    osr[:,:,1] = x ⊗ i2
    osr[:,:,2] = x ⊗ i2

    soq = i2 ⊗ x + x ⊗ i2 + x ⊗ x + x ⊗ x
    sor = zeros(dtype,(4,4,2))
    sor[:,:,1] = i2 ⊗ z
    sor[:,:,2] = i2 ⊗ z

    ooq = i2 ⊗ x + x ⊗ i2 + z ⊗ x + z ⊗ x
    ool = zeros(dtype,(4,4,2))
    ool[:,:,1] = i2 ⊗ z; ool[:,:,2] = i2 ⊗ z
    oor = zeros(dtype,(4,4,2))
    oor[:,:,1] = x ⊗ i2; oor[:,:,2] = x ⊗ i2
    oop = zeros(dtype,(4,4,2,2))

    sos = -(i2 ⊗ osq + x ⊗ i2 ⊗ i2 + x ⊗ osr[:,:,1] + x ⊗ osr[:,:,2])

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
    @test isapprox(s*o*s, sos)
end
