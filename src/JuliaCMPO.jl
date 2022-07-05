module JuliaCMPO
__precompile__()

using Dates, Parameters, Reexport  
using CUDA; CUDA.allowscalar(false)
using Random; Random.seed!()
using HDF5, DelimitedFiles, Printf
using LogExpFunctions
using OMEinsum, LinearAlgebra, KrylovKit
using Zygote, Optim, ChainRules
@reexport using FiniteTLanczos

@reexport import Base: kron, *, eltype, size, length, getindex, iterate
@reexport import Base: ==, ≈,  transpose, adjoint, cat
@reexport import LinearAlgebra: ishermitian, norm, normalize, diag, diagm
@reexport import KrylovKit.eigsolve
@reexport import FiniteTLanczos.eigensolver
@reexport import FiniteTLanczos: partitian,
                                 free_energy,
                                 energy,
                                 specific_heat,
                                 entropy,
                                 thermal_average, 
                                 correlation_2time

# structs/CMPSandCMPO
export AbstractCTensor, AbstractCMPS, AbstractCMPO,
       CMPS, CMPO, CuCMPS, CuCMPO,
       CMPS_generate, CMPO_generate,
       CTensor, CuCTensor,
       bond_dimension,
       virtual_dimension
#structs/CMPSMatrix
export CMPSMatrix

#utilities
export ⊗, diagm

# logtrexp.jl
export logtrexp

# solver
export cpu_solver, gpu_solver
# OptimFunctions
export veclength, optim_functions


# CMPSOperations
export log_overlap,
       logfidelity, fidelity, 
       project,
       diagQ
# CMPSInitiate
export init_cmps, 
       interpolate_isometry, 
       MeraUpdateStep, MeraUpdateTrace,
       MeraUpdateResult, MeraUpdateOptions,
       adaptive_mera_update,
       CompressResult, 
       compress_cmps

# PhysicalModels
export PhysModel, Ising_CMPO, generalUt, expand_cmpo
export TFIsing, TFIsing_2D_helical,
       XYmodel, XYmodel_2D_helical, 
       XXZmodel,  XXZmodel_2D_helical

       
# SaveAndLoad
export saveCMPS, readCMPS

# ThermaldynamicQuanties
export make_operator

# Correlations
export Lehmann_spectral_function, Lehmann_A, 
       Lehmann_structure_factor, Lehmann_S

# evaluate
export evaluate, evaluate_py, 
       hermitian_evaluate, 
       non_hermitian_evaluate


include("./structs/CMPSandCMPO.jl")
include("./structs/CMPSMatrix.jl")

include("./utilities/otimes.jl")
include("./utilities/diagm.jl")

include("eigensolver.jl")
include("logtrexp.jl")

include("solver.jl")
include("OptimFunctions.jl")

include("CTensorProducts.jl")
include("CMPSOperations.jl")
include("CMPSInitiate.jl")

include("PhysicalModels.jl")
include("SaveAndLoad.jl")

include("ThermaldynamicQuanties.jl")
include("Correlations.jl")

include("evaluate.jl")
include("./rrule/logtrexp.jl")


end # module
