using JuliaCMPO, Test
processor = CPU
solver = solver_function(processor)
@show solver

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

@testset "gradient/logtrexp.jl" begin
    include("./gradient/logtrexp.jl")
end