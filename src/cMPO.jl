module cMPO
__precompile__()

using LinearAlgebra, GenericLinearAlgebra
using Zygote, FiniteDifferences, ChainRulesCore
using Optim, FluxOptTools
using Random; Random.seed!()
using StatsFuns, SpecialFunctions, HCubature
using OMEinsum

import Base: *, isequal, transpose, adjoint, cat

# utilities
export pauli,
       delta,
       Masubara_freq,
       symmetrize, symeigen, 
       logtrexp

#SaveAndLoad
export saveCMPS, readCMPS

# structs
export CMPS, CMPO, PhysModel

# operations
export toarray, tovector, tocmps,
       #normalize, 
       log_overlap,
       transpose, adjoint, ishermitian,
       project,
       diagQ

# multiplications
export ⊗

# cMPSCompress
export init_cmps, cmps_compress

# PhysicalModels
export Ising_CMPO, generalUt, expand
export TFIsing, XYmodel, XXZmodel, #HeisenbergModel,
       TFIsing_2D_helical,
       XYmodel_2D_helical,
       XXZmodel_2D_helical

# ThermaldynamicQuanties
export make_operator,
       thermal_average,
       free_energy,
       energy,
       specific_heat,
       entropy

# Correlations
export correlation_2time,
       LehmannGFKernel, Masubara_freq_GF,
       Lehmann_spectral_function, Lehmann_A, 
       Lehmann_structure_factor, Lehmann_S

# evaluate
export evaluate,
       hermitian_evaluate,
       non_hermitian_evaluate


include("utilities.jl")
include("SaveAndLoad.jl")
include("structs.jl")
include("multiplications.jl")
include("operations.jl")
include("cMPSCompress.jl")
include("PhysicalModels.jl")
include("ThermaldynamicQuanties.jl")
include("Correlations.jl")
include("evaluate.jl")

#include("rrule.jl")
#inclide("PowerProjection.jl")

end # module
