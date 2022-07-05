using JuliaCMPO, Test
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
        D1 = 4;  D2 = 2
        T = ComplexF64
        a = [randn(D1,D1),randn(T,D1,D1),randn(D1,D1,D2), randn(T,D1,D1,D2),randn(D1,D1,D2,D2),randn(T,D1,D1,D2,D2)]
        for i = 1:length(a), j = i+1 : min(2*(div(i+1,2)+1), length(a))
            @test Array(solver(⊗, a[i], a[j])) ≈ otimes(a[i],a[j])
        end
    end
end
