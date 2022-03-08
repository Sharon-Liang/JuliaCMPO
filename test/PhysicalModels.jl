using Test, cMPO
using LinearAlgebra

@testset "Ising_CMPO" begin
    J = 1.0
    pz = pauli(:z); px = pauli(:x)
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

@testset "cat two cmpos" begin
    pz = pauli(:z); px = pauli(:x)
    for w1 in [1,2,3], w2 in [1,2,3]
        o1 = Ising_CMPO(1.0, pz, pz, w1)
        o2 = Ising_CMPO(1.0, px, px, w2)

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

@testset "expand a cmpo" begin
    pz = pauli(:z); px = pauli(:x)
    for w1 in [1, 3]
        o1 = Ising_CMPO(1.0, pz, pz, w1)
        o2 = Ising_CMPO(0.5, pz, pz, w1)

        O1 = expand(o1)
        O2 = cat(o2, transpose(o2))

        @test isequal(O1, O2)
    end
end

@testset "Ut' * T * Ut = transpose(T)" begin
    wid = 3
    for m in [TFIsing(1.0,1.0), XYmodel(), XYmodel_2D_helical(1), XXZmodel_2D_helical(2.0, wid)]
        T = m.Tmatrix
        Ut = m.Ut
        @test isequal(Ut' * T * Ut, transpose(T))
    end
end