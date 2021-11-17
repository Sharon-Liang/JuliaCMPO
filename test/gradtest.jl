using cMPO, Test
using Random; Random.seed!()
using Zygote, ForwardDiff, FiniteDifferences

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

function gradcheck(f, xs...)
    grad_zygote = Zygote.gradient(f, xs...)
    grad_finite_difference = ngradient(f, xs...)
    return all(isapprox.(grad_zygote, grad_finite_difference; rtol = 1e-5, atol = 1e-5))
end

gradtest(f, xs::AbstractArray...) = gradcheck((xs...) -> sum(sin.(f(xs...))), xs...)
gradtest(f, dims...) = gradtest(f, rand.(Float64, dims)...)


@testset "⊗" begin
    D1 = 3;  D2 = 4
    T = Float64
    a = [randn(D1,D1),randn(T,D1,D1),randn(D1,D1,D2),randn(T,D1,D1,D2),randn(D1,D1,D2,D2),randn(T,D1,D1,D2,D2)]
    for i = 1:length(a), j = i+1 : min(2*(div(i+1,2)+1), length(a))
        @test gradtest(⊗, a[i], a[j])
    end
end



@testset "Gradient test: TFIsing model" begin
    β = 20.0; χ = 2
    g = rand(1)[1]
    w = TFIsing(1.,g)
    ψ = init_cmps(χ) |> toarray
    fenergy = ψ->free_energy(ψ, w, β)
    @test gradcheck(fenergy, ψ)
end

"""
@testset "Gradient test: Heisenberg model" begin
    β = 1.0; χ = 2
    W = HeisenbergModel()
    ψ = init_cmps(χ,D=3) |> toarray
    fenergy = ψ->free_energy(ψ, W, β)
    @test gradcheck(fenergy, ψ)
end
"""