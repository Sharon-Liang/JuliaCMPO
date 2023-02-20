using Test, JuliaCMPO
using JLD, DelimitedFiles

processor = CPU
obsv_functions = [free_energy, energy, entropy]


#=
### *TFIsingChain* : *Variational Principle*
=#
βlist = [1.0, 1.5]
Tₘ = model(TFIsingChain(), 1.0)

#exact values
Oe = zeros(2, 4)
Oe[1, :] = [1.00, -1.4152076398462621, -1.1179418373401746, 0.2972658025060875]
Oe[2, :] = [1.50, -1.3338698801984272, -1.2073855733208403, 0.1897264603163804]


@testset "TFIsingChain_variation" begin
    bondD = 8
    result_folder = "./test/TFIsingChain_variation"
    
    variation_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder)
    Oc = readdlm("./test/TFIsingChain_variation/obsvs.txt", skipstart = 1)
    @test Oc ≈ Oe rtol = 1.e-5
end


@testset "TFIsingChain_power" begin
    bondD = 8
    result_folder = "./test/TFIsingChain_power"
    
    power_evaluate(Tₘ, bondD, βlist; processor, obsv_functions, result_folder, max_pow_step = 50)
    Oc = readdlm("./test/TFIsingChain_variation/obsvs.txt", skipstart = 1)
    @test Oc ≈ Oe rtol = 1.e-3

    #test start from cMPS
    result_folder = "./test/TFIsingChain_power_init"
    init_step = 40
    init = Vector(undef, 2)
    init[1] = load("test/TFIsingChain_power/cmps_ckpt_beta_1.000000.jld")[string(init_step)]
    init[2] = nothing
    
    power_evaluate(Tₘ, bondD, βlist, init; processor, obsv_functions, result_folder, max_pow_step = 50)
    Oc = readdlm("./test/TFIsingChain_variation/obsvs.txt", skipstart = 1)
    @test Oc ≈ Oe rtol = 1.e-3
end