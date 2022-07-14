using JuliaCMPO, Test
using CUDA; CUDA.allowscalar(false)
mydevice = GPU
solver = solver_function(mydevice)

@testset "otimes.jl" begin
    include("otimes.jl")
end

@testset "CTensorProducts.jl" begin
    include("CTensorProducts.jl")
end

@testset "CMPSMatrix.jl" begin
    include("CMPSMatrix.jl")
end