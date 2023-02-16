using Test, JuliaCMPO
using LinearAlgebra, Parameters
using Random; Random.seed!()
using FiniteDifferences, Zygote

ngradient(f, x) = grad(central_fdm(5, 1), f, x)
zgradient = Zygote.gradient

solver = cpu_solver
#=
### *logtrexp* 
=#
@testset "logtrexp-matrix" begin
    M = rand(10, 10)
    t = rand()
    f = x -> logtrexp(t, x)

    ngrad = ngradient(f, M)[1]
    zgrad = zgradient(f, M)[1]
    @test ngrad ≈ zgrad
end


#=
### *logfidelity* 
=#
@testset "logfidelity" begin
    β = rand()
    χ = 4
    for D = 1:2
        ψ₀ = init_cmps(χ+2, D) |> solver
        ψ  = init_cmps(χ, D) |> diagQ |> solver

        Qd = Vector(diag(ψ.Q))
        R  = Array(ψ.R)

        pars = Zygote.Params([Qd, R])
        loss() = -logfidelity(solver(CMPS(diagm(Qd), R)), ψ₀, β)

        p₀, f, g! = optim_functions(loss, pars)

        ngrad = ngradient(f, p₀)[1]
        zgrad = similar(p₀); g!(zgrad, p₀)
        
        @test zgrad ≈ ngrad
    end
end


#=
### *project function* 
=#
@testset "project function" begin 
    β = rand()
    χ = 6
    for D = 1:2
        u = rand(χ, χ-2) |>  solver
        ψ₀ = init_cmps(χ, D) |> solver

        loss(u) = logfidelity(project(ψ₀, u), ψ₀, β, false)

        ngrad = ngradient(loss, u)[1]
        zgrad = zgradient(loss, u)[1]

        @test zgrad ≈ ngrad
    end
end