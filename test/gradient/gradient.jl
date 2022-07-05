using JuliaCMPO, Test
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
    include("./logtrexp.jl")
end


@testset "logfidelity" begin
    include("logfidelity.jl")
end

@testset "project function" begin
    include("./project.jl")
end

@testset "free energy" begin
    include("./free_energy.jl")
end













