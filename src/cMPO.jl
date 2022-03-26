module cMPO
__precompile__()

using Random; Random.seed!()
using HDF5, DelimitedFiles, Printf
using StatsFuns
using OMEinsum, LinearAlgebra #, GenericLinearAlgebra
using Zygote, Optim
using FiniteDifferences

import Base: *, isequal, transpose, adjoint, cat
import LinearAlgebra: ishermitian, norm, normalize

# utilities
export pauli,
       delta,
       Masubara_freq,
       symmetrize, symeigen, 
       logtrexp

# OptimFunctions
export veclength, optim_functions,
       optim_functions_py

# structs
export CMPS, CMPO, PhysModel, 
       MeraUpdateStep, MeraUpdateTrace,
       MeraUpdateResult, MeraUpdateOptions,
       CompressResult

# SaveAndLoad
export saveCMPS, readCMPS

# cMPSOperations
export toarray, tovector, tocmps,
       norm, normalize, 
       log_overlap,
       transpose, adjoint, ishermitian,
       project,
       diagQ

# multiplications
export ⊗

# cMPSInitiate
export init_cmps, 
       logfidelity, fidelity, 
       interpolate_isometry, adaptive_mera_update,
       compress_cmps
       #init_cmps_py,
       #compress_cmps_py

# PhysicalModels
export Ising_CMPO, generalUt, expand_cmpo
export TFIsing, 
       XYmodel, 
       XXZmodel, 
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
export evaluate, evaluate_py, 
       hermitian_evaluate, #hermitian_evaluate_py,
       non_hermitian_evaluate#, non_hermitian_evaluate_py


include("utilities.jl")
include("OptimFunctions.jl")
include("structs.jl")
include("SaveAndLoad.jl")
include("multiplications.jl")
include("cMPSOperations.jl")
include("cMPSInitiate.jl")
include("PhysicalModels.jl")
include("ThermaldynamicQuanties.jl")
include("Correlations.jl")
include("evaluate.jl")

#include("rrule.jl")

end # module
