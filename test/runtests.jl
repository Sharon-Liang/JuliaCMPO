using JuliaCMPO, Test

filelist = ["util.jl", "math.jl", "core.jl", "grads.jl"]


processor = CPU
for file in filelist
    include(file)
end

@testset "models.jl" begin
    include("models.jl")
end


