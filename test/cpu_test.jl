using JuliaCMPO, Test
solver = cpu_solver

@testset "otimes.jl" begin
    include("otimes.jl")
end

@testset "CTensorProducts.jl" begin
    include("CTensorProducts.jl")
end

@testset "CMPSMatrix.jl" begin
    include("CMPSMatrix.jl")
end