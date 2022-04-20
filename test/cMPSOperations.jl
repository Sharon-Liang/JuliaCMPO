using Test, cMPO
using LinearAlgebra
using Random; Random.seed!()


for solver in [cpu_solver, gpu_solver] 
    @testset "$(solver)" begin
        @testset "normalize CMPS" begin
            β = rand(1:5)
            ψ = init_cmps(rand(2:5), rand(1:2))

            ψ = solver(x->normalize(x, β), ψ)
            @test solver(x->norm(x, β),ψ) ≈ 1.
        end

        @testset "Ut' * T * Ut = transpose(T)" begin
            wid = 3
            for m in [TFIsing(1.0,1.0), XYmodel(), XYmodel_2D_helical(1, expand=true), XXZmodel_2D_helical(2.0, wid, expand=true)]
                Tm = solver(x->x, m.Tmatrix)
                Ut = solver(x->x, m.Ut)
                @test  Ut' * Tm * Ut == transpose(Tm)
            end
        end
    end
end


@testset "adjoint and ishermitian" begin
    for m in [TFIsing(1.0,1.0), XYmodel(), XXZmodel(1.0)]
        T = m.Tmatrix
        @test ==(T, T')
        @test ishermitian(T)
    end

    wid = 3
    for m in [TFIsing_2D_helical(1.0,1.0, wid), XYmodel_2D_helical(wid), XXZmodel_2D_helical(1.0, wid)]
        T = m.Tmatrix
        @test ==(T, T') == false
        @test ishermitian(T) == false
    end
end

@testset "cat two CMPOs" begin
    pz = pauli(PZ); px = pauli(PX)
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

@testset "expand a CMPO" begin
    pz = pauli(PZ); px = pauli(PX)
    for w1 in [1, 3]
        o1 = Ising_CMPO(1.0, pz, pz, w1)
        o2 = Ising_CMPO(0.5, pz, pz, w1)

        O1 = expand_cmpo(o1)
        O2 = cat(o2, transpose(o2))

        @test ==(O1, O2)
    end
end

@testset "transpose(T^2) = transpose(T) * transpose(T)" begin
    #one should compare transfer matrix instead of cMPO local tensor
    β = 2.0
    for m in [TFIsing(1.0,1.0), XXZmodel_2D_helical(1.0, 2)]
        ψ = init_cmps(2, m.vir_dim)
        T1 =  transpose(m.Tmatrix * m.Tmatrix)
        T2 =  transpose(m.Tmatrix) * transpose(m.Tmatrix)
        @test free_energy(ψ, T1, β) ≈ free_energy(ψ, T2, β)
    end
end