using JuliaCMPO, Test
using CUDA; CUDA.allowscalar(false)
processor = GPU
solver = solver_function(processor)

@show solver

@testset "otimes.jl" begin
    include("otimes.jl")
end

@testset "CTensorProducts.jl" begin
    include("CTensorProducts.jl")
end

@testset "CMPSMatrix.jl" begin
    include("CMPSMatrix.jl")
end

@testset "logtrexp.jl" begin
    include("logtrexp.jl")
end