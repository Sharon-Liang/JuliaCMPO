using Test, JuliaCMPO
using LinearAlgebra, Parameters
using Random; Random.seed!()
using FiniteDifferences, Zygote

ngradient(f, x) = grad(central_fdm(5, 1), f, x)
zgradient = Zygote.gradient

@show solver = solver_function(processor)

#=
### *logtrexp* 
=#
@testset "logtrexp-matrix" begin
    M = rand(10, 10)
    t = rand()
    f = x -> logtrexp(t, x)

    ngrad = ngradient(f, M)[1]
    zgrad = zgradient(f, solver(M))[1]
    @test ngrad ≈ zgrad |> Array
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
        u = rand(χ, χ-2) 
        ψ₀ = init_cmps(χ, D) |> solver

        loss(u) = logfidelity(project(ψ₀, solver(u)), ψ₀, β, false)

        ngrad = ngradient(loss, u)[1]
        zgrad = zgradient(loss, u)[1]

        @test zgrad ≈ ngrad
    end
end



#=
### *free energy* 
=#
@testset "Free energy: TFIsing Chain" begin
    β = 2.0
    χ = 4
    m = model(TFIsingChain(), 1.0) |> solver
    
    ψ = init_cmps(χ) |> diagQ

    Qd = convert(Vector, diag(ψ.Q))
    R = convert(Array, ψ.R)
    pars = Zygote.Params([Qd, R])
    loss() = free_energy(solver(CMPS(diagm(Qd), R)), m, β)
    
    p₀, f, g! = optim_functions(loss, pars)

    ngrad = ngradient(f, p₀)[1]
    zgrad = similar(p₀); g!(zgrad, p₀)
        
    @test zgrad ≈ ngrad
end



#=
### *shift_spectrum* 
=#
@testset "shift_spectrum" begin
    β = 2.0
    χ = 4
    Tₘ = model(TFIsingChain(), 1.0) |> solver

    to_shift = 1.e-3

    ψ₀ = init_cmps(χ+2) |> solver
    ψ  = init_cmps(χ) |> diagQ |> solver

    Qd = Vector(diag(ψ.Q))
    R  = Array(ψ.R)

    pars = Zygote.Params([Qd, R])
    function loss()
        ψ = solver(CMPS(diagm(Qd), R))
        Ol₁ = log_overlap(ψ, Tₘ * ψ₀, β) |> exp
        Ol₂ = log_overlap(ψ, ψ₀, β) |> exp
        N₀ = 0.5 * log_overlap(ψ, ψ, β)
        return -log(Ol₁ + to_shift * Ol₂) + N₀
    end

    p₀, f, g! = optim_functions(loss, pars)

    ngrad = ngradient(f, p₀)[1]
    zgrad = similar(p₀); g!(zgrad, p₀)
        
    @test zgrad ≈ ngrad
end
