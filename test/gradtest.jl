#https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl

function ngradient(f, xs::AbstractArray...) #finite difference
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
    grad_zygote = gradient(f, xs...)
    grad_finite_difference = ngradient(f, xs...)
    return all(isapprox.(grad_zygote, grad_finite_difference; rtol = 1e-5, atol = 1e-5))
end

gradtest(f, xs::AbstractArray...) = gradcheck((xs...) -> sum(sin.(f(xs...))), xs...)
#gradtest(f, xs::AbstractArray...) = gradcheck((xs...) -> sum( abs.(sin.(f(xs...)))), xs...)
gradtest(f, dims...) = gradtest(f, rand.(Float64, dims)...)


Dtype = [Float64] #Fail for ComplexF32
for dtype in Dtype
    str = @sprintf "Input data type = %s" dtype
    @testset "$str" begin
        @testset "Test Gradient of ⊗" begin
            a1 = rand(dtype,3,3); a2 = rand(3,3)
            b1 = rand(dtype,3,3,3); b2 = rand(3,3)
            c1 = rand(dtype,3,3,3,3); c2 = rand(3,3,3,3)
            @test gradtest(⊗, a1, a2)
            @test gradtest(⊗, a1, b1)
            @test gradtest(⊗, b1, a1)
            @test gradtest(⊗, b1, b2)
            @test gradtest(⊗, b1, c1)
            @test gradtest(⊗, c1, b1)
            @test gradtest(⊗, c1, c2)
        end
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