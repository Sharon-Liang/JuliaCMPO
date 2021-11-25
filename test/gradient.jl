using cMPO, Test
using Random; Random.seed!()
using Zygote, ForwardDiff, FiniteDiff

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


function gradcheck_fdiff(f, xs)
    grad_fdiff = ForwardDiff.gradient(f, xs)
    grad_finite_difference = ngradient(f, xs)[1]
    return all(isapprox.(grad_fdiff, grad_finite_difference; rtol = 1e-5, atol = 1e-5))
end
gradtest_fdiff(f, xs::AbstractArray) = gradcheck_fdiff(xs -> sum(sin.(f(xs))), xs)

@testset "⊗: zygote" begin
    D1 = 3;  D2 = 4
    T = Float64
    a = [randn(D1,D1),randn(T,D1,D1),randn(D1,D1,D2),randn(T,D1,D1,D2),randn(D1,D1,D2,D2),randn(T,D1,D1,D2,D2)]
    for i = 1:length(a)
        @test gradtest(x->x⊗x, a[i])
        for j = i+1: min(2*(div(i+1,2)+1), length(a))
            @test gradtest(⊗, a[i], a[j])
        end
    end
end

function test_prod(sl::AbstractArray, sr::AbstractArray)
    sl = sl |> tocmps
    sr = sr |> tocmps
    sl * sr
end

function test_prod(sl::AbstractArray, o::CMPO; sr::AbstractArray)
    sl = sl |> tocmps
    sr = sr |> tocmps
    sl * o * sr
end

function test_prod(o::CMPO, sr::AbstractArray; sl::AbstractArray)
    sl = sl |> tocmps
    sr = sr |> tocmps
    sl * o * sr
end

function test_prod(sl::AbstractArray, o::CMPO, sr::AbstractArray)
    sl = sl |> tocmps
    sr = sr |> tocmps
    sl * o * sr
end

@testset "*: zygote" begin
    χ = 4
    sl = init_cmps(χ) |> toarray
    sr = init_cmps(χ) |> toarray
    o = TFIsing(1.0, 1.0)
    @test gradtest(test_prod, sl, sr)
    @test gradtest(x -> test_prod(x,o; sr = sr), sl)
    @test gradtest(x -> test_prod(o,x; sl = sl), sr)
    @test gradtest((x,y) -> test_prod(x,o,y), sl, sr)
end

@testset "*: forwarddiff" begin
    χ = 4
    sl = init_cmps(χ) |> toarray
    sr = init_cmps(χ) |> toarray
    o = TFIsing(1.0, 1.0)
    @test gradtest_fdiff(x -> test_prod(x, sr) , sl)
    @test gradtest_fdiff(x -> test_prod(sl, x) , sr)
    @test gradtest_fdiff(x -> test_prod(x,o; sr = sr), sl)
    @test gradtest_fdiff(x -> test_prod(o,x; sl = sl), sr)
    @test gradtest(x -> test_prod(x,o,sr), sl)
    @test gradtest(x -> test_prod(sl,o,x), sr)
end

@testset "logtrexp" begin
    D = 8
    A = rand(D, D) |> symmetrize
    @test gradtest(logtrexp, A)  
    @test gradtest_fdiff(logtrexp, A)  
end

@testset "free energy: TFIsing model" begin
    β = 20.0; χ = 2
    g = rand(1)[1]
    w = TFIsing(1.,g)
    ψ = init_cmps(χ) |> toarray
    @test gradcheck(x->free_energy(x, w, β), ψ)
    @test gradcheck_fdiff(x->free_energy(x, w, β), ψ)
end

@testset "Gradient test: XY model" begin
    β = 1.0; χ = 2
    w = XYmodel()
    ψ = init_cmps(χ,D=2) |> toarray
    @test gradcheck(ψ->free_energy(ψ, w, β), ψ)
    @test gradcheck_fdiff(x->free_energy(x, w, β), ψ)
end

@testset "Gradient test: AFM Heisenberg model" begin
    β = 1.0; χ = 2
    w = HeisenbergModel()
    ψ = init_cmps(χ,D=3) |> toarray
    @test gradcheck(ψ->free_energy(ψ, w, β), ψ)
    @test gradcheck_fdiff(x->free_energy(x, w, β), ψ)
end

@testset "Gradient test: XXZ model" begin
    β = 1.0; χ = 2
    g = rand(1)[1]
    w = XXZmodel(g)
    g == 0 ?  dR = 2 : dR = 3
    ψ = init_cmps(χ,D=dR) |> toarray
    @test gradcheck(ψ->free_energy(ψ, w, β), ψ)
    @test gradcheck_fdiff(x->free_energy(x, w, β), ψ)
end