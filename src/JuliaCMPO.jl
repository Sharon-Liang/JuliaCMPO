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
@reexport import LinearAlgebra: ishermitian, norm, normalize, diag
@reexport import KrylovKit.eigsolve
@reexport import FiniteTLanczos: eigensolver, symmetrize, 
                                 Processor, CPU, GPU, 
                                 solver_function,
                                 cpu_solver, gpu_solver,
                                 FTLMOptions, TraceEstimator


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

#utilities
export ⊗


# OptimFunctions.jl
export veclength, optim_functions
# SaveLoad.jl
export saveCMPS, readCMPS

# logtrexp.jl
export logtrexp
# CMPSOperations
export log_overlap,
       logfidelity, fidelity, 
       project,
       diagQ

#MeraUpdate.jl
export interpolate_isometry, adaptive_mera_update
#CMPSCompress.jl 
export compress_cmps
# CMPSInitiate
export init_cmps

#PhysicalQuantities/PhysicalModels.jl
export PhysModel, Ising_CMPO, generalUt, expand_cmpo
export TFIsing, TFIsing_2D_helical,
       XYmodel, XYmodel_2D_helical, 
       XXZmodel,  XXZmodel_2D_helical

#PhysicalQuantities/Thermaldynamic.jl
export make_operator

#structs/EvaluateOptions.jl
export EstimatorType, EvaluateOptions

#evaluate.jl
#export evaluate,
#       hermitian_evaluate, 
#       non_hermitian_evaluate


include("./structs/CMPSandCMPO.jl")
include("./structs/CMPSMatrix.jl")
include("./structs/MeraUpdate.jl")
include("./structs/CMPSCompress.jl")


include("./utilities/otimes.jl")

include("solver.jl")
include("OptimFunctions.jl")
include("SaveLoad.jl")

include("CTensorProducts.jl")

include("eigensolver.jl")
include("logtrexp.jl")
include("CMPSOperations.jl")

include("MeraUpdate.jl")
include("CMPSCompress.jl")
include("CMPSInitiate.jl")

include("./PhysicalQuantities/PhysicalModels.jl")
include("./PhysicalQuantities/Thermaldynamic.jl")
include("./PhysicalQuantities/Correlations.jl")

include("./structs/EvaluateOptions.jl")
include("evaluate.jl")

include("./rrule/ConstructorAdjoint.jl")
include("./rrule/accum.jl")
include("./rrule/logtrexp.jl")


end # module
