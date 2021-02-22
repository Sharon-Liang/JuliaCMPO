using Test
using Zygote
using Optim
using LinearAlgebra
using Random ; Random.seed!()
using cMPO
"""
Gradient test
"""
@testset "Zygote gradient" begin
    β = 1.
    χ = 2
    W = TFIsing(1.,1.)
    ψ = init_cmps(χ)

    # Numerical Difference
    ϵ = 1.e-5
    g = zeros(χ, χ, 2)

    for i = 1:length(ψ.Q)
        ψ.Q[i] -= ϵ
        f1 = FreeEnergy(ψ, W, β)

        ψ.Q[i] += 2ϵ
        f2 = FreeEnergy(ψ, W, β)

        g[i] = (f2 - f1)/ 2ϵ
        ψ.Q[i] -= ϵ
    end

    for i = 1:length(ψ.R)
        ψ.R[i] -= ϵ
        f1 = FreeEnergy(ψ, W, β)

        ψ.R[i] += 2ϵ
        f2 = FreeEnergy(ψ, W, β)

        g[i+χ^2] = (f2 - f1)/ 2ϵ
        ψ.R[i] -= ϵ

    end
        g2 = Zygote.gradient(ψ->FreeEnergy(ψ, W, β), ψ)[1]
        g2 = cmps(g2.Q, g2.R) |> toarray
        for i = 1: length(g2)
           @test isapprox(g[i], g2[i], rtol=1e-2)
        end
end
