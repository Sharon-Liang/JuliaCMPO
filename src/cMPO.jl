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

export pauli, ⊗

export symmetrize, trexp, value, logtrexp
export grad_func, grad_num
export cmps, cmpo
export toarray, init_cmps, ovlp

export TFIsing
export free_energy
export thermal_average, correlation_2time
export susceptibility, imag_susceptibility, structure_factor
export energy_density, ave_sx, critical_zz_cor, critical_zz_sus

include("Setup.jl")
include("PhysicalObservables.jl")
include("utilities.jl")

include("exact.jl")

#@reexport using .Setup
#@reexport using .PhysicalObservables

end # module
