using cMPO, Test

@testset "setup" begin
    include("setup.jl")
end

@testset "multiplication" begin
    include("multiplication.jl")
end

@testset "gradient" begin
    include("gradient.jl")
end

@testset "hessian" begin
    include("hessian.jl")
end
