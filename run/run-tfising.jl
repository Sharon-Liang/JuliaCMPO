using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

println("Start!")
#include("TFIsing-change-gamma.jl")
include("Sw_compare.jl")

