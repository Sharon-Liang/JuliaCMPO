module cMPO
__precompile__()

using LinearAlgebra, GenericLinearAlgebra
using Zygote, FiniteDifferences
using Optim
using Random; Random.seed!()
using StatsFuns, SpecialFunctions, HCubature
using OMEinsum

import Base: *, isequal
import LinearAlgebra: normalize

export pauli, 
       delta, 
       Masubara_freq, 
       ⊗
export symmetrize, 
       TrExp, trexp, value, logtrexp
export gradient_function, hessian_function
export CMPS, CMPO
export toarray, tovector, tocmps
export init_cmps, ovlp

export make_operator
export TFIsing, XYmodel, XXZmodel, HeisenbergModel
export thermal_average, partitian, partitian!
export free_energy, energy, specific_heat, entropy
export correlation_2time,
       Masubara_freq_GF


include("rrule.jl")
include("utilities.jl")
include("setup.jl")
include("PhysicalModels.jl")
include("PhysicalObservables.jl")


end # module
