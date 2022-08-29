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

@reexport import FiniteTLanczos: eigensolver, symmetrize, 
                                 Processor, CPU, GPU, 
                                 solver_function,
                                 cpu_solver, gpu_solver,
                                 FTLMOptions, TraceEstimator


# structs/CMPSandCMPO.jl
export AbstractCTensor, AbstractCMPS, AbstractCMPO,
       CMPS, CMPO, CuCMPS, CuCMPO,
       cmps_generate, cmpo_generate,
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
#export update_processor

#utilities
export ⊗


# OptimFunctions.jl
export veclength, optim_functions
# SaveLoad.jl
export saveCMPS, readCMPS

# logtrexp.jl
export logtrexp
# CMPSOperations
export tomatrix,
       log_overlap,
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
export ising_cmpo
export TFIsing, TFIsing_2D_helical,
       XYmodel, XYmodel_2D_helical, 
       XXZmodel, XXZmodel_2D_helical

#PhysicalQuantities/Thermaldynamic.jl
export make_operator

#structs/EvaluateOptions.jl
export EstimatorType, EvaluateOptions


include("./structs/CMPSandCMPO.jl")
include("./structs/CMPSMatrix.jl")
include("./structs/MeraUpdate.jl")
include("./structs/CMPSCompress.jl")


include("./utilities/otimes.jl")

include("solver.jl")
include("optim_functions.jl")
include("save_load.jl")

include("ctensor_products.jl")

include("eigensolver.jl")
include("logtrexp.jl")
include("cmps_operations.jl")

include("mera_update.jl")
include("cmps_compress.jl")
include("cmps_initiate.jl")

include("./physical_quantities/physical_models.jl")
include("./physical_quantities/thermaldynamic.jl")
#include("./physical_quantities/Correlations.jl")

include("./structs/EvaluateOptions.jl")
include("evaluate.jl")

include("./rrule/logtrexp.jl")


end # module
