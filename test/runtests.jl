using JuliaCMPO, Test

@testset "cpu_test.jl" begin
    include("cpu_test.jl")
end

@testset "gpu_test.jl" begin
    include("gpu_test.jl")
end




