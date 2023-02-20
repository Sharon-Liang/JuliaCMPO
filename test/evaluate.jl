using Test, JuliaCMPO
using JLD, DelimitedFiles
using LinearAlgebra

processor = CPU
obsv_functions = [free_energy, energy, entropy]



### *TFIsingChain*
βlist = [1.0, 1.5]
Tₘ = model(TFIsingChain(), 1.0)

#exact values
Oe = zeros(2, 4)
Oe[1, :] = [1.00, -1.4152076398462621, -1.1179418373401746, 0.2972658025060875]
Oe[2, :] = [1.50, -1.3338698801984272, -1.2073855733208403, 0.1897264603163804]


#=
### *TFIsingChain* : *Variational Principle*
=#
#=
@testset "TFIsingChain_variation" begin
    bondD = 8
    result_folder = "./test/TFIsingChain_variation"
    
    variation_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder)
    Oc = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Oc ≈ Oe rtol = 1.e-5
end
=#



#=
### *TFIsingChain* : *Power Method*
=#
@testset "TFIsingChain_power" begin
    bondD = 8
    result_folder = "./test/TFIsingChain_power"
    
    power_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder, max_pow_step = 50)
    Oc = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Oc ≈ Oe rtol = 1.e-3

    #test start from cMPS
    init_step = 40
    init = Vector(undef, 2)
    init[1] = load(result_folder*"/cmps_ckpt_beta_1.000000.jld")[string(init_step)]
    init[2] = nothing
    
    new_result_folder = "./test/TFIsingChain_power_init"
    power_evaluate(Tₘ, bondD, βlist, init; processor, obsv_functions, result_folder = new_result_folder, max_pow_step = 50)
    Oc = readdlm(new_result_folder*"/obsvs.txt", skipstart = 1)
    @test Oc ≈ Oe rtol = 1.e-3
end


#=
### *HeisenbergChain* : *Power method*
=#
βlist = [1.0, 2.0]
obsv_functions = [free_energy]
Tₘ = model(XXZChain(), 1.0)

#dtrg result with bond dimension 256
Oe = zeros(2, 2)
Oe[:, 1] = [1.000000, -0.7953882215614290]
Oe[:, 2] = [2.000000, -0.5385576248164411] 
@testset "HeisenbergChain_power" begin
    bondD = 8
    result_folder = "./test/HeisenbergChain_power"
    
    power_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder, max_pow_step = 50)
    Oc = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Oc ≈ Oe rtol = 1.e-3

    to_group = 2
    result_folder = "./test/HeisenbergChain_power_group"
    power_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder, max_pow_step = 50, to_group)
    Og = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Og ≈ Oe rtol = 1.e-3
    @test norm(Og - Oe) < norm(Oc - Oe)
end
