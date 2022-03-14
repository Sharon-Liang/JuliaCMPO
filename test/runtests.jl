using cMPO, Test

@testset "utilities.jl" begin
    include("utilities.jl")
end

@testset "multiplication.jl" begin
    include("multiplications.jl")
end

@testset "operations.jl" begin
    include("operations.jl")
end

@testset "cmpsInitiate.jl" begin
    include("cmpsInitiate.jl")
end

#@testset "gradient" begin
#    include("gradient.jl")
#end

#@testset "hessian" begin
#    include("hessian.jl")
#end
