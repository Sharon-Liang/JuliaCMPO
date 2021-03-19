using LinearAlgebra
using StatsFuns
using DelimitedFiles
using JLD, HDF5

include("exact.jl")

# exact-sx
β = 20
len = 51
g = [i for i in range(0,2,length = len)]

open("./data/Exact_sx.txt","w") do io
    for j in g
        sx = ave_sx(1., j, β)
        writedlm(io,[j sx])
    end
end
