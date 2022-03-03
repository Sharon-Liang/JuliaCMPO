module cMPO
__precompile__()

using LinearAlgebra, GenericLinearAlgebra
using Zygote, FiniteDifferences, ChainRulesCore
using Optim
using Random; Random.seed!()
using StatsFuns, SpecialFunctions, HCubature
using OMEinsum

import Base: *, isequal, transpose, adjoint
import LinearAlgebra: normalize

export pauli,
       delta,
       Masubara_freq,
       symmetrize, symeigen, 
       logtrexp,
       gradient_function, hessian_function

export CMPS, CMPO,
       toarray, tovector, tocmps,
       normalize, 
       log_overlap,
       transpose,
       project,
       diagQ

export ⊗

export init_cmps

export PhysModel, TFIsing,
       XYmodel, XXZmodel, HeisenbergModel

export make_operator,
       thermal_average,
       free_energy,
       energy,
       specific_heat,
       entropy

export correlation_2time,
       LehmannGFKernel, Masubara_freq_GF,
       Lehmann_spectral_function, Lehmann_A, 
       Lehmann_structure_factor, Lehmann_S


include("utilities.jl")
include("setup.jl")
include("cMPSInitiate.jl")
include("rrule.jl")
include("multiplications.jl")
inclide("PowerProjection.jl")
include("PhysicalModels.jl")
include("ThermaldynamicQuanties.jl")
include("Correlations.jl")



end # module
