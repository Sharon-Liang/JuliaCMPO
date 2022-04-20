using cMPO, Test
using LinearAlgebra
using Random; Random.seed!()

function otimes(A::Matrix, B::Matrix)
    kron(A,B)
end

function otimes(A::Matrix, B::Array{T,3} where T)
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(B)[3]
    eltype(A) == eltype(B) ? dtype = eltype(A) : dtype = ComplexF64
    res = zeros(dtype, χ1*χ2, χ1*χ2, D)
    for d = 1:D
        res[:,:,d] = otimes(A, B[:,:,d])
    end
    return res
end

function otimes(A::Array{T,3} where T, B::Matrix)
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(A)[3]
    eltype(A) == eltype(B) ? dtype = eltype(A) : dtype = ComplexF64
    res = zeros(dtype, χ1*χ2, χ1*χ2, D)
    for d = 1:D
        res[:,:,d] = otimes(A[:,:,d], B)
    end
    return res
end

function otimes(A::Array{T,3} where T, B::Array{T,3} where T)
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(A)[3]
    eltype(A) == eltype(B) ? dtype = eltype(A) : dtype = ComplexF64
    res = zeros(dtype, χ1*χ2, χ1*χ2)
    for d = 1:D
        res += otimes(A[:,:,d], B[:,:,d])
    end
    return res
end

function otimes(A::Array{T,4} where T, B::Array{T,3} where T)
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(B)[3]
    eltype(A) == eltype(B) ? dtype = eltype(A) : dtype = ComplexF64c
    res = zeros(dtype, χ1*χ2, χ1*χ2,D)
    for d = 1:D
        res[:,:,d] = otimes(A[:,:,d,:], B)
    end
    return res
end

function otimes(A::Array{T,3} where T, B::Array{T,4} where T)
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(A)[3]
    eltype(A) == eltype(B) ? dtype = eltype(A) : dtype = ComplexF64
    res = zeros(dtype,χ1*χ2, χ1*χ2, D)
    for d = 1:D
        res[:,:,d] = otimes(A, B[:,:,:,d])
    end
    return res
end

function otimes(A::Array{T,4} where T, B::Array{T,4} where T)
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(A)[3]
    eltype(A) == eltype(B) ? dtype = eltype(A) : dtype = ComplexF64
    res = zeros(dtype, χ1*χ2, χ1*χ2, D, D)
    for d1 = 1:D, d2 = 1:D
        res[:,:,d1, d2] = otimes(A[:,:,d1,:], B[:,:,:,d2])
    end
    return res
end


for solver in [cpu_solver, gpu_solver] 
    @testset "$(solver)" begin
        @testset "⊗" begin
            D1 = 4;  D2 = 2
            T = ComplexF64
            a = [randn(D1,D1),randn(T,D1,D1),randn(D1,D1,D2), randn(T,D1,D1,D2),randn(D1,D1,D2,D2),randn(T,D1,D1,D2,D2)]
            for i = 1:length(a), j = i+1 : min(2*(div(i+1,2)+1), length(a))
                @test Array(solver(⊗, a[i], a[j])) ≈ otimes(a[i],a[j])
            end
        end

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
                
                s = solver(CMPS_generate, x, z)
                o = solver(CMPO_generate, x, x, z, p)

                @test  s*s ≈ solver(A->A, ss)         
                @test (o*s).Q ≈ solver(A->A, osq)
                @test (o*s).R ≈ solver(A->A, osr)
                @test (s*o).Q ≈ solver(A->A, soq) 
                @test (s*o).R ≈ solver(A->A, sor) 
                @test (o*o).Q ≈ solver(A->A, ooq)  
                @test (o*o).R ≈ solver(A->A, oor) 
                @test (o*o).L ≈ solver(A->A, ool) 
                @test (o*o).P ≈ solver(A->A, oop) 
                @test s*o*s ≈ solver(A->A, sos)
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
            
                s = solver(CMPS_generate, x, R)
                o = solver(CMPO_generate, x, R, L, P)

                @test  s*s ≈ solver(A->A, ss)         
                @test (o*s).Q ≈ solver(A->A, osq)
                @test (o*s).R ≈ solver(A->A, osr)
                @test (s*o).Q ≈ solver(A->A, soq) 
                @test (s*o).R ≈ solver(A->A, sor) 
                @test (o*o).Q ≈ solver(A->A, ooq)  
                @test (o*o).R ≈ solver(A->A, oor) 
                @test (o*o).L ≈ solver(A->A, ool) 
                @test (o*o).P ≈ solver(A->A, oop) 
                @test s*o*s ≈ solver(A->A, sos)
            end
        end
    end
end
