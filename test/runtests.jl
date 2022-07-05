using JuliaCMPO, Test
using CUDA; CUDA.allowscalar(false)

@testset "utilities.jl" begin
    include("utilities.jl")
end

@testset "otimes.jl" begin
    include("otimes.jl")
end

@testset "CTensorProducts.jl" begin
    include("CTensorProducts.jl")
end

@testset "CMPSMatrix.jl" begin
    include("CMPSMatrix.jl")
end

@testset "CMPSOperations.jl" begin
    include("CMPSOperations.jl")
end

@testset "CMPSInitiate.jl" begin
    include("CMPSInitiate.jl")
end

@testset "PhysicalModels.jl" begin
    include("PhysicalModels.jl")
end

@testset "gradient" begin
    include("./gradient/gradient.jl")
end




