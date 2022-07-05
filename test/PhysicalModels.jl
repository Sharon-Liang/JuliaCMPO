using Test, JuliaCMPO
using LinearAlgebra

@testset "Ising_CMPO" begin
    J = 1.0
    pz = pauli(PZ); px = pauli(PX)
    # H = -J px pz
    @testset "width = 1" begin
        o = Ising_CMPO(J, px, pz)
        @test o.Q == zeros(2,2)
        @test o.R == pz
        @test o.L == px
        @test o.P == zeros(2,2)
    end

    @testset "width > 1" begin
        wid = 5
        o = Ising_CMPO(J, px, pz, wid)
        @test o.Q == zeros(2,2)
        r = zeros(2,2,wid)
        r[:,:,1] = pz
        @test o.R == r
        l = zeros(2,2,wid)
        l[:,:,1] = px; l[:,:,wid] = px
        @test o.L == l
        p = zeros(2,2,wid,wid)
        for i = 2:wid
            p[:,:, i, i-1] = Matrix{Float64}(I,2,2)
        end
        @test o.P == p
    end
end

