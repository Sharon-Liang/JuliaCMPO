using cMPO, Test
using Random; Random.seed!()
using Zygote
using LinearAlgebra

function ngradient(f, xs::AbstractArray...)
    #https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl
    grads = zero.(xs)
    for (x, Δ) in zip(xs, grads), i in 1:length(x)
      δ = sqrt(eps())
      tmp = x[i]
      x[i] = tmp - δ/2
      y1 = f(xs...)
      x[i] = tmp + δ/2
      y2 = f(xs...)
      x[i] = tmp
      Δ[i] = (y2-y1)/δ
    end
    return grads
end

@testset "logtrexp" begin
    D = rand(2:6)
    A = rand(D, D)

    grad_ndiff = ngradient(logtrexp, A)
    for device in [:cpu, :gpu]
        grad_zygote = Zygote.gradient(M->logtrexp(M,device=device), A)
        @test all(isapprox.(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5))
    end
end

@testset "logfidelity" begin
    β = 2.0
    χ = 6
    vir_dim = rand(1:2)
    ψ0 = init_cmps(χ+2, vir_dim)
    ψ = init_cmps(χ, vir_dim)

    for device in [:cpu, :gpu]
        loss() = -logfidelity(CMPS(diag(ψ.Q)|> diagm, ψ.R), ψ0, β, device= device)
        p0, f, g! = optim_functions(loss, Params([ψ.Q, ψ.R]))

        grad_ndiff = ngradient(f, p0)[1]
        grad_zygote = similar(p0); g!(grad_zygote, p0)
        @test all(isapprox.(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5))
    end
end

@testset "free energy: TFIsing model" begin
    β = 2.0; χ = 2
    g = rand(1)[1]
    m = TFIsing(1.,g)
    ψ = init_cmps(χ) |> diagQ
    for device in [:cpu, :gpu]
        loss() = free_energy(CMPS(diag(ψ.Q)|> diagm, ψ.R), m.Tmatrix, β)
        p0, f, g! = optim_functions(loss, Params([ψ.Q, ψ.R]))
        
        grad_ndiff = ngradient(f, p0)[1]
        grad_zygote = similar(p0); g!(grad_zygote, p0)
        @test all(isapprox.(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5))
    end
end












