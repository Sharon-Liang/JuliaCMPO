using Test, JuliaCMPO
using LinearAlgebra
using Random; Random.seed!()

@show solver = solver_function(processor)

#=
### *Multiplications* : *Otimes* 
=#
println("Define otimes functions")
otimes(A::Matrix, B::Matrix) = kron(A,B)

function otimes(A::Matrix{T}, B::Array{T,3}) where T
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(B)[3]
    res = zeros(T, χ1*χ2, χ1*χ2, D)
    for d = 1:D
        res[:,:,d] = otimes(A, B[:,:,d])
    end
    return res
end

function otimes(A::Array{T,3}, B::Matrix{T}) where T
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(A)[3]
    res = zeros(T, χ1*χ2, χ1*χ2, D)
    for d = 1:D
        res[:,:,d] = otimes(A[:,:,d], B)
    end
    return res
end

function otimes(A::Array{T,3}, B::Array{T,3}) where T
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(A)[3]
    res = zeros(T, χ1*χ2, χ1*χ2)
    for d = 1:D
        res += otimes(A[:,:,d], B[:,:,d])
    end
    return res
end

function otimes(A::Array{T,4}, B::Array{T,3}) where T
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(B)[3]
    res = zeros(T, χ1*χ2, χ1*χ2,D)
    for d = 1:D
        res[:,:,d] = otimes(A[:,:,d,:], B)
    end
    return res
end

function otimes(A::Array{T,3}, B::Array{T,4}) where T
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(A)[3]
    res = zeros(T, χ1*χ2, χ1*χ2, D)
    for d = 1:D
        res[:,:,d] = otimes(A, B[:,:,:,d])
    end
    return res
end

function otimes(A::Array{T,4}, B::Array{T,4}) where T
    χ1 = size(A)[1]; χ2 = size(B)[1]; D = size(A)[3]
    res = zeros(T, χ1*χ2, χ1*χ2, D, D)
    for d1 = 1:D, d2 = 1:D
        res[:,:,d1, d2] = otimes(A[:,:,d1,:], B[:,:,:,d2])
    end
    return res
end


println("Start testing...")
@testset "otimes" begin
    D1 = 4;  D2 = 2
    a = [rand(D1,D1), rand(D1,D1,D2), rand(D1,D1,D2,D2)]
    
    for i in eachindex(a), j in eachindex(a)
        if abs(lastindex(size(a[i]))  - lastindex(size(a[j]))) < 2
            val₁ = solver(⊗, a[i], a[j])
            val₂ = otimes(a[i],a[j]) |> solver
            @test val₁ ≈ val₂ 
        end
    end
end


#=
### *Multiplications* : *Products between CMPS and CMPO* 
=#

@testset "multiplications of cmps and cmpo: D-1 = 1" begin
    D = 4

    x = symmetrize(rand(D,D))
    z = symmetrize(rand(D,D))
    p = zeros(D,D)
                
    i2 = oneunit(x)
    ss = -(i2 ⊗ x + x ⊗ i2 + z ⊗ z)

    osq = i2 ⊗ x + x ⊗ i2 + z ⊗ z
    osr = x ⊗ i2
                    
    soq = i2 ⊗ x + x ⊗ i2 + z ⊗ x
    sor = i2 ⊗ z

    ooq = i2 ⊗ x + x ⊗ i2 + z ⊗ x
    ool = i2 ⊗ z
    oor = x ⊗ i2
    oop = zeros(D^2, D^2)

    sos = -(i2 ⊗ osq + x ⊗ i2 ⊗ i2 + z ⊗ osr)
                
    s = solver(CMPS(x, z))
    o = solver(CMPO(x, x, z, p))

    @test  s*s ≈ solver(ss)         
    @test  o*s ≈ solver(CMPS(osq, osr))
    @test  s*o ≈ solver(CMPS(soq, sor)) 
    @test  o*o ≈ solver(CMPO(ooq, oor, ool, oop))
    @test  s*o*s ≈ solver(sos)  
end

@testset "multiplications of cmps and cmpo: D-1 > 1" begin
    D = 4

    x = symmetrize(rand(D,D))
    z = symmetrize(rand(D,D))
    p = zeros(D,D,2,2)
                
    i2 = oneunit(x)
    ss = -(i2 ⊗ x + x ⊗ i2 + x ⊗ x + x ⊗ x)

    osq = i2 ⊗ x + x ⊗ i2 + z ⊗ x + z ⊗ x
    osr = cat(x ⊗ i2, x ⊗ i2, dims=3)

    soq = i2 ⊗ x + x ⊗ i2 + x ⊗ x + x ⊗ x
    sor = cat(i2 ⊗ z, i2 ⊗ z, dims=3)

    ooq = i2 ⊗ x + x ⊗ i2 + z ⊗ x + z ⊗ x
    ool = cat(i2 ⊗ z, i2 ⊗ z, dims=3)
    oor = cat(x ⊗ i2, x ⊗ i2, dims=3)
    oop = zeros(D^2, D^2, 2, 2)

    sos = -(i2 ⊗ osq + x ⊗ i2 ⊗ i2 + x ⊗ osr[:,:,1] + x ⊗ osr[:,:,2])

    R = cat(x, x, dims=3)
    L = cat(z, z, dims=3)
    P = zeros(D,D,2,2)

    s = solver(CMPS(x, R))
    o = solver(CMPO(x, R, L, P))

    @test  s*s ≈ solver(ss)         
    @test  o*s ≈ solver(CMPS(osq, osr))
    @test  s*o ≈ solver(CMPS(soq, sor)) 
    @test  o*o ≈ solver(CMPO(ooq, oor, ool, oop))  
    @test  s*o*s ≈ solver(sos) 
end


#=
### *Properties* : *CMPS* 
=#
@testset "normalize a CMPS" begin
    β = rand()
    for D = 1:2
        ψ = init_cmps(10, D) |> solver
        ψ1 = normalize(ψ, β)
        @test norm(ψ1, β) ≈ 1.
    end
end


