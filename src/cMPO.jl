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

export pauli, delta, Masubara_freq, ⊗

export symmetrize, trexp, value, logtrexp
export gradient_function, hessian_function
export CMPS, CMPO
export toarray, tovector, tocmps
export init_cmps, ovlp

export make_operator
export TFIsing, XYmodel, XXZmodel, HeisenbergModel
export thermal_average, partitian, partitian!
export free_energy, energy, specific_heat, entropy
export correlation_2time
export check_anomalous_term, Masubara_freq_GF, Masubara_freq_GFdivOmega
export spectral_density, susceptibility, structure_factor
export energy_density, ave_sx

include("rrule.jl")
include("utilities.jl")
include("setup.jl")
include("PhysicalModels.jl")
include("PhysicalObservables.jl")

include("exact.jl")


end # module
