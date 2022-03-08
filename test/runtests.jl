using cMPO, Test

@testset "utilities.jl" begin
    include("utilities.jl")
end

@testset "multiplication" begin
    include("multiplications.jl")
end

@testset "gradient" begin
    include("gradient.jl")
end

@testset "hessian" begin
    include("hessian.jl")
end
