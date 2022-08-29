using LinearAlgebra
using Random; Random.seed!()

@testset "multiplications of cmps and cmpo: D-1 = 1" begin
    T = ComplexF64
    D = 4
    a = [randn(D,D), randn(T,D,D)]
    for i = 1:length(a), j=1:length(a)
        x = symmetrize(a[i])
        z = symmetrize(a[j])
        p = zeros(T,D,D)
                
        i2 = Matrix{T}(I,D,D)
        ss = -(i2 ⊗ x + x ⊗ i2 + z ⊗ z)

        osq = i2 ⊗ x + x ⊗ i2 + z ⊗ z
        osr = x ⊗ i2
                    
        soq = i2 ⊗ x + x ⊗ i2 + z ⊗ x
        sor = i2 ⊗ z

        ooq = i2 ⊗ x + x ⊗ i2 + z ⊗ x
        ool = i2 ⊗ z
        oor = x ⊗ i2
        oop = zeros(T,D^2,D^2)

        sos = -(i2 ⊗ osq + x ⊗ i2 ⊗ i2 + z ⊗ osr)
                
        s = solver(cmps_generate, x, z)
        o = solver(cmpo_generate, x, x, z, p)

        @test  tomatrix(s*s) ≈ solver(ss)         
        @test (o*s).Q ≈ solver(osq)
        @test (o*s).R ≈ solver(osr)
        @test (s*o).Q ≈ solver(soq) 
        @test (s*o).R ≈ solver(sor) 
        @test (o*o).Q ≈ solver(ooq)  
        @test (o*o).R ≈ solver(oor) 
        @test (o*o).L ≈ solver(ool) 
        @test (o*o).P ≈ solver(oop) 
        @test tomatrix(s*o*s) ≈ solver(sos)
    end
end

@testset "multiplications of cmps and cmpo: D-1 > 1" begin
    T = ComplexF64
    D = 4
    a = [randn(D,D), randn(T,D,D)]
    for i = 1:length(a), j=1:length(a)
        x = symmetrize(a[i])
        z = symmetrize(a[j])
        p = zeros(T,D,D,2,2)
                
        i2 = Matrix{T}(I,D,D)
        ss = -(i2 ⊗ x + x ⊗ i2 + x ⊗ x + x ⊗ x)

        osq = i2 ⊗ x + x ⊗ i2 + z ⊗ x + z ⊗ x
        osr = zeros(T, D^2, D^2, 2)
        osr[:,:,1] = x ⊗ i2
        osr[:,:,2] = x ⊗ i2

        soq = i2 ⊗ x + x ⊗ i2 + x ⊗ x + x ⊗ x
        sor = zeros(T, D^2, D^2, 2)
        sor[:,:,1] = i2 ⊗ z
        sor[:,:,2] = i2 ⊗ z

        ooq = i2 ⊗ x + x ⊗ i2 + z ⊗ x + z ⊗ x
        ool = zeros(T, D^2, D^2, 2)
        ool[:,:,1] = i2 ⊗ z; ool[:,:,2] = i2 ⊗ z
        oor = zeros(T, D^2, D^2, 2)
        oor[:,:,1] = x ⊗ i2; oor[:,:,2] = x ⊗ i2
        oop = zeros(T, D^2, D^2, 2, 2)

        sos = -(i2 ⊗ osq + x ⊗ i2 ⊗ i2 + x ⊗ osr[:,:,1] + x ⊗ osr[:,:,2])

        R = zeros(T,D,D,2); R[:,:,1] = x ; R[:,:,2] = x
        L = zeros(T,D,D,2); L[:,:,1] = z ; L[:,:,2] = z
        P = zeros(T,D,D,2,2)
            
        s = solver(cmps_generate, x, R)
        o = solver(cmpo_generate, x, R, L, P)

        @test  tomatrix(s*s) ≈ solver(ss)         
        @test (o*s).Q ≈ solver(osq)
        @test (o*s).R ≈ solver(osr)
        @test (s*o).Q ≈ solver(soq) 
        @test (s*o).R ≈ solver(sor) 
        @test (o*o).Q ≈ solver(ooq)  
        @test (o*o).R ≈ solver(oor) 
        @test (o*o).L ≈ solver(ool) 
        @test (o*o).P ≈ solver(oop) 
        @test tomatrix(s*o*s) ≈ solver(sos)
    end
end

