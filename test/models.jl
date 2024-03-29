using Test, JuliaCMPO
using LinearAlgebra

import JuliaCMPO: _ising_2D_block

@testset "Ising_CMPO" begin
    J = 1.0
    pz = pauli(PZ); px = pauli(PX)
    # H = -J px pz
    @testset "width = 1" begin
        o = _ising_2D_block(J, px, pz)
        @test o.Q == zeros(2,2)
        @test o.R == pz
        @test o.L == px
        @test o.P == zeros(2,2)
    end

    @testset "width > 1" begin
        wid = 5
        o = _ising_2D_block(J, px, pz, wid)
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

@testset "cat two CMPOs" begin
    pz = pauli(PZ); px = pauli(PX)
    for w1 in [1,2,3], w2 in [1,2,3]
        o1 = _ising_2D_block(1.0, pz, pz, w1)
        o2 = _ising_2D_block(1.0, px, px, w2)

        o = cat(o1, o2)
        @test o.Q == zeros(2,2)

        r = zeros(2,2,w1+w2)
        r[:,:,1] = pz
        r[:,:,w1+1] = px
        @test o.R == r

        l = zeros(2,2,w1+w2)
        l[:,:,1] = pz; l[:,:,w1] = pz
        l[:,:,w1+1] = px; l[:,:,w1+w2] = px
        @test o.L == l

        p = zeros(2,2,w1+w2, w1+w2)
        for i = 2:w1
            p[:,:, i, i-1] = Matrix{Float64}(I,2,2)
        end
        for i = w1+2 : w1+w2
            p[:,:, i, i-1] = Matrix{Float64}(I,2,2)
        end
        @test o.P == p 
    end   
end

@testset "adjoint and ishermitian" begin
    for M in [TFIsingChain(1.0, 1.0), XXChain(), XXZChain(1.0)]  
        T = model(M)
        @test ==(T, T')
        @test ishermitian(T)
    end

    W = 3
    for M in [TFIsingSquareHelical(1.0,1.0,W), XXSquareHelical(W), XXZSquareHelical(1.0,W)]   
        T = model(M)
        @test ==(T, T') == false
        @test ishermitian(T) == false
    end
end


@testset "transpose(T²) = transpose(T) * transpose(T)" begin
    #one should compare transfer matrix instead of JuliaCMPO local tensor
    β = 2.0
    for M in [TFIsingChain(1.0,1.0), XXZSquareHelical(1.0, 2)]
        T = model(M)
        ψ = init_cmps(2, T, β)
        T1 =  transpose(T * T)
        T2 =  transpose(T) * transpose(T)
        @test free_energy(ψ, T1, β) ≈ free_energy(ψ, T2, β)
    end
end

