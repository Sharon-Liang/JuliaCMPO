using Test
using Zygote
using Optim
using LinearAlgebra
using Random ; Random.seed!()
using cMPO


@testset "setup_test.jl" begin
    include("setup_test.jl")
end

@testset "gradtest.jl" begin
    include("gradtest.jl")
end

