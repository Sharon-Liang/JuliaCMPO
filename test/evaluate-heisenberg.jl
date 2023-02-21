using Test, JuliaCMPO
using JLD, DelimitedFiles
using LinearAlgebra

processor = CPU
obsv_functions = [free_energy]


#=
### *HeisenbergChain* : *Power method*
=#
βlist = [1.0, 2.0]
Tₘ = model(XXZChain(), 1.0)

#dtrg result with bond dimension 256
Oe = zeros(2, 2)
Oe[1, :] = [1.000000, -0.7953882215614290]
Oe[2, :] = [2.000000, -0.5385576248164411] 


#=
### *HeisenbergChain* : *Power method*
=#
#=
@testset "HeisenbergChain_power" begin
    bondD = 8
    result_folder = "./HeisenbergChain_power"
    
    power_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder, max_pow_step = 50)
    Oc = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Oc ≈ Oe rtol = 1.e-3
end 


@testset "HeisenbergChain_power_group" begin
    bondD = 8
    result_folder = "./HeisenbergChain_power_group"

    to_group = 2
    
    power_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder, max_pow_step = 50, to_group)
    Og = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Og ≈ Oe rtol = 1.e-3
end
=#

@testset "HeisenbergChain_power_shift" begin
    bondD = 8
    result_folder = "./HeisenbergChain_power_shift"
    
    to_shift = 1.e-3
    
    power_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder, max_pow_step = 50, to_shift)
    Os = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Os ≈ Oe rtol = 1.e-3
end
