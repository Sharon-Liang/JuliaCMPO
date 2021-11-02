module cMPO
__precompile__()

using LinearAlgebra
using Zygote, FiniteDifferences
using Optim
using Random; Random.seed!()
using StatsFuns, SpecialFunctions
using OMEinsum

import Base: *, isequal
import LinearAlgebra: normalize

export pauli, delta, Masubara_freq, ⊗

export symmetrize, trexp, value, logtrexp
export gradient_function, hessian_function
export cmps, cmpo
export toarray, tovector, tocmps
export init_cmps, ovlp

export make_operator
export TFIsing, XYmodel, HeisenbergModel
export thermal_average, partitian, partitian!
export free_energy, energy, specific_heat, entropy
export correlation_2time
export check_anomalous_term, Masubara_freq_GF, Masubara_freq_GFdivOmega
export spectral_density, susceptibility, structure_factor
export energy_density, ave_sx


include("utilities.jl")
include("setup.jl")
include("PhysicalModels.jl")
include("PhysicalObservables.jl")

include("exact.jl")


end # module
