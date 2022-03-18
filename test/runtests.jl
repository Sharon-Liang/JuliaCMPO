using cMPO, Test

@testset "utilities.jl" begin
    include("utilities.jl")
end

@testset "multiplication.jl" begin
    include("multiplications.jl")
end

@testset "cMPSOperations.jl" begin
    include("cMPSOperations.jl")
end

@testset "cMPSInitiate.jl" begin
    include("cMPSInitiate.jl")
end

@testset "gradient" begin
    include("gradient.jl")
end

#@testset "hessian" begin
#    include("hessian.jl")
#end
