using JuliaCMPO, Test
processor = CPU
solver = solver_function(processor)

@testset "PhysicalModels.jl" begin
    include("PhysicalModels.jl")
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

@testset "logtrexp.jl" begin
    include("./logtrexp.jl")
end