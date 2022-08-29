using JuliaCMPO, Test
using CUDA; CUDA.allowscalar(false)
processor = GPU
solver = solver_function(processor)
solver = solver_function(processor)

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

@testset "CMPSOperations.jl" begin
    include("CMPSOperations.jl")
end

#gradient test
@testset "gradient/logtrexp.jl" begin
    include("./gradient/logtrexp.jl")
end

@testset "gradient/logfidelity.jl" begin
    include("./gradient/logfidelity.jl")
end

@testset "gradient/project.jl" begin
    include("./gradient/project.jl")
end

@testset "gradient/free_energy.jl" begin
    include("gradient/free_energy.jl")
end

@testset "CMPSInitiate.jl" begin
    include("CMPSInitiate.jl")
end