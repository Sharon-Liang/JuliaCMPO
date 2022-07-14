using JuliaCMPO, Test
using LinearAlgebra, Random; Random.seed!()
using Zygote, FiniteDifferences

ngradient(f, x) = grad(central_fdm(5, 1), f, x)
zgradient = Zygote.gradient

@testset "logtrexp-matrix" begin
    include("./logtrexp.jl")
end


#@testset "logfidelity" begin
#    include("logfidelity.jl")
#end

#@testset "project function" begin
#    include("./project.jl")
#end

#@testset "free energy" begin
#    include("./free_energy.jl")
#end













