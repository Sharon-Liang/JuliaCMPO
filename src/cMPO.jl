module cMPO
__precompile__()

using LinearAlgebra, GenericLinearAlgebra
using Zygote, FiniteDifferences, ChainRulesCore
using Optim
using Random; Random.seed!()
using StatsFuns, SpecialFunctions, HCubature
using OMEinsum

import Base: *, isequal, transpose, adjoint, cat
import LinearAlgebra: normalize

# utilities
export pauli,
       delta,
       Masubara_freq,
       symmetrize, symeigen, 
       logtrexp,
       gradfunc_gen, hessfunc_gen

# structs
export CMPS, CMPO, PhysModel

# operations
export toarray, tovector, tocmps,
       #normalize, 
       log_overlap,
       transpose, adjoint, 
       project,
       diagQ

# multiplications
export ⊗

# cMPSInitiate
export init_cmps

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


include("utilities.jl")
include("structs.jl")
include("operations.jl")
include("cMPSInitiate.jl")
include("rrule.jl")
include("multiplications.jl")
#inclide("PowerProjection.jl")
include("PhysicalModels.jl")
include("ThermaldynamicQuanties.jl")
include("Correlations.jl")



end # module
