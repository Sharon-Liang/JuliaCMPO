using JuliaCMPO, Test
using CUDA; CUDA.allowscalar(false)
device = CPU
solver = solver_function(device)

@testset "otimes.jl" begin
    include("otimes.jl")
end

@testset "CTensorProducts.jl" begin
    include("CTensorProducts.jl")
end

@testset "CMPSMatrix.jl" begin
    include("CMPSMatrix.jl")
end