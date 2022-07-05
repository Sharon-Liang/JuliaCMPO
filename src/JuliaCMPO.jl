module JuliaCMPO
__precompile__()

using Dates, Parameters, Reexport  
using Random; Random.seed!()
using HDF5, DelimitedFiles, Printf
using LogExpFunctions
using OMEinsum, LinearAlgebra, KrylovKit
using Zygote, Optim, ChainRules

using CUDA; CUDA.allowscalar(false)

import Base: kron, *, eltype, size, length, getindex, iterate
import Base: ==, ≈,  transpose, adjoint, cat
import LinearAlgebra: ishermitian, norm, normalize, diag, diagm
import KrylovKit.eigsolve

using FiniteTLanczos
import FiniteTLanczos: eigensolver,
                       partitian,
                       free_energy,
                       energy,
                       specific_heat,
                       entropy,
                       thermal_average

#eigensolver.jl
export eigensolver

# structs/CMPSandCMPO
export AbstractCTensor, AbstractCMPS, AbstractCMPO,
       CMPS, CMPO, CuCMPS, CuCMPO,
       CMPS_generate, CMPO_generate,
       CTensor, CuCTensor,
       bond_dimension,
       virtual_dimension

#structs/CMPSMatrix
export CMPSMatrix


# utilities
export PauliMatrixName, PX, PY, iPY, PZ, PPlus, PMinus,
       pauli,
       delta,
       OperatorType, Fermi, Bose,
       Masubara_freq,
       symmetrize, symeigen, 
       logtrexp,
       diagm

# solver
export cpu_solver, gpu_solver

# OptimFunctions
export veclength, optim_functions


# multiplications
export ⊗

# TraceEstimator/cMPSMatrix.jl


# TraceEstimator/trexpsolver.jl
export eigsolve

# cMPSOperations
export log_overlap,
       norm, normalize, 
       logfidelity, fidelity, 
       ishermitian,
       project,
       diagQ

# PhysicalModels
export PhysModel, Ising_CMPO, generalUt, expand_cmpo
export TFIsing, TFIsing_2D_helical,
       XYmodel, XYmodel_2D_helical, 
       XXZmodel,  XXZmodel_2D_helical

# cMPSInitiate
export init_cmps, 
       interpolate_isometry, 
       MeraUpdateStep, MeraUpdateTrace,
       MeraUpdateResult, MeraUpdateOptions,
       adaptive_mera_update,
       CompressResult, 
       compress_cmps


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


# SaveAndLoad
export saveCMPS, readCMPS

# evaluate
export evaluate, evaluate_py, 
       hermitian_evaluate, 
       non_hermitian_evaluate

include("./structs/CMPSandCMPO.jl")
include("./structs/CMPSMatrix.jl")

include("utilities.jl")

include("eigensolver.jl")


include("solver.jl")
include("OptimFunctions.jl")

include("multiplications.jl")
include("cMPSOperations.jl")




include("PhysicalModels.jl")
include("cMPSInitiate.jl")

include("ThermaldynamicQuanties.jl")
include("Correlations.jl")

include("SaveAndLoad.jl")
include("evaluate.jl")
include("rrule.jl")


end # module
