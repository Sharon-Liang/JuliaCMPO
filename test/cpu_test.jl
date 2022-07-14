using JuliaCMPO, Test
mydevice = CPU
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

@testset "logtrexp.jl" begin
    include("./logtrexp.jl")
end