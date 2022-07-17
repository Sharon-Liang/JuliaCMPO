using JuliaCMPO, Test
using FiniteDifferences, Zygote

ngradient(f, x) = grad(central_fdm(5, 1), f, x)
zgradient = Zygote.gradient

@testset "cpu_test.jl" begin
    include("cpu_test.jl")
end

@testset "gpu_test.jl" begin
    include("gpu_test.jl")
end




