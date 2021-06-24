module cMPO
__precompile__()

#using Reexport:@reexport
using LinearAlgebra
using Zygote
using Optim
using Random; Random.seed!()
using StatsFuns, SpecialFunctions

import Base: *, isequal
import LinearAlgebra: normalize

export pauli, delta, ⊗

export symmetrize, trexp, value, logtrexp
export grad_func, grad_num
export cmps, cmpo
export toarray, init_cmps, ovlp, tocmps

export TFIsing, XYmodel, HeisenbergModel
export thermal_average, partitian, partitian!
export free_energy, energy, specific_heat, entropy
export correlation_2time
export susceptibility, imag_susceptibility, structure_factor
export energy_density, ave_sx, critical_zz_cor, critical_zz_sus
export critical_zz_chi

include("Setup.jl")
include("PhysicalObservables.jl")
include("utilities.jl")

include("exact.jl")

#@reexport using .Setup
#@reexport using .PhysicalObservables

end # module
