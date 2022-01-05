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
       symmetrize, 
       TrExp, trexp, value, logtrexp,
       gradient_function, hessian_function

export CMPS, CMPO,
       toarray, tovector, tocmps,
       init_cmps, 
       normalize, ovlp

export ⊗

export TFIsing, 
       XYmodel, XXZmodel, HeisenbergModel

export make_operator,
       thermal_average, 
       partitian,
       free_energy, 
       energy, 
       specific_heat, 
       entropy,
       correlation_2time,
       Masubara_freq_GF


include("rrule.jl")
include("utilities.jl")
include("setup.jl")
include("multiplication.jl")
include("PhysicalModels.jl")
include("PhysicalObservables.jl")


end # module
