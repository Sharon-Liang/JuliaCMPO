module JuliaCMPO
__precompile__()

using Dates, Parameters, Reexport  
using CUDA; CUDA.allowscalar(false)
using Random; Random.seed!()
using HDF5, DelimitedFiles, Printf
using LogExpFunctions
using OMEinsum, LinearAlgebra, KrylovKit
using Zygote, Optim, ChainRules
using PhysModels, FiniteTLanczos

@reexport import Base: kron, *, eltype, size, length, getindex, iterate
@reexport import Base: ==, ≈,  transpose, adjoint, cat
@reexport import LinearAlgebra: ishermitian, norm, normalize, diag, diagm
@reexport import KrylovKit.eigsolve
@reexport import FiniteTLanczos: eigensolver, symmetrize


# structs/CMPSandCMPO.jl
export AbstractCTensor, AbstractCMPS, AbstractCMPO,
       CMPS, CMPO, CuCMPS, CuCMPO,
       CMPS_generate, CMPO_generate,
       CTensor, CuCTensor,
       bond_dimension,
       virtual_bond_dimension

#structs/CMPSMatrix.jl
export CMPSMatrix

#structs/MeraUpdate.jl
export MeraUpdateOptions,
       MeraUpdateStep,
       MeraUpdateTrace,
       MeraUpdateResult

#structs/CMPSCompress.jl
export CompressOptions, CompressResult

#structs/EvaluateOptions.jl
export EstimatorType, EvaluateOptions

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
       adaptive_mera_update,
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

# evaluate
export evaluate, evaluate_py, 
       hermitian_evaluate, 
       non_hermitian_evaluate

export PhysModel 

include("./structs/CMPSandCMPO.jl")
include("./structs/CMPSMatrix.jl")
include("./structs/MeraUpdate.jl")
include("./structs/CMPSCompress.jl")


include("./utilities/otimes.jl")
include("./utilities/diagm.jl")

include("solver.jl")
include("OptimFunctions.jl")
include("SaveLoad.jl")

include("eigensolver.jl")
include("logtrexp.jl")

include("CTensorProducts.jl")
include("CMPSOperations.jl")
include("MeraUpdate.jl")
include("CMPSCompress.jl")
include("CMPSInitiate.jl")

include("./PhysicalQuantities/PhysicalModels.jl")
include("./PhysicalQuantities/Thermaldynamic.jl")
include("./PhysicalQuantities/Correlations.jl")

include("./structs/EvaluateOptions.jl")
include("evaluate.jl")


include("./rrule/logtrexp.jl")


end # module
