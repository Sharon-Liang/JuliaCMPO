using Test
using Zygote
using Optim
using LinearAlgebra
using Random ; Random.seed!()
using cMPO

@testset "setup_test.jl" begin
    include("setup_test.jl")
end

@testset "Gradient test" begin
    β = 1.0; χ = 2
    W = TFIsing(1.,1.)
    ψ = init_cmps(χ) |> toarray
    f_eng = ψ->free_energy(ψ, W, β)
    g_num = grad_num(f_eng, ψ)
    g_auto = gradient(f_eng, ψ)[1]
    @test isapprox(g_auto[:,:,1], g_num[:,:,1], rtol=1e-5)
    @test isapprox(g_auto[:,:,2], g_num[:,:,2], rtol=1e-5)
end
