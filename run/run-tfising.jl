using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf
using TimerOutputs

include("TFIsing-change-gamma.jl")
include("hessian-data.jl")

