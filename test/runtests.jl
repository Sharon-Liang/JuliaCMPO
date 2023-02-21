using JuliaCMPO, Test
using FiniteDifferences, Zygote

#=
filelist = ["util.jl", "math.jl", "core.jl", "grads.jl"]


processor = CPU
for file in filelist
    @testset "$(processor): $(file)" begin
        include(file)
    end
end

@testset "models.jl" begin
    include("models.jl")
end
=#



filelist = ["core.jl"]
processor = GPU
for file in filelist
    @testset "$(processor): $(file)" begin
        include(file)
    end
end


#include("evaluate-heisenberg.jl")