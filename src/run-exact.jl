using LinearAlgebra
using StatsFuns
using DelimitedFiles
using JLD, HDF5

include("exact.jl")

# exact free energy
Γ = 2.0
len = 40
beta = [i for i in range(1,20,length = len)]

open("./data/0323/fexact-gamma-2.0.txt","w") do io
    for β in beta
        f = free_energy(1., Γ, β)
        writedlm(io,[β f])
    end
end


# exact sx
open("./data/0323/sx-exact-gamma-2.0.txt","w") do io
    for β in beta
        sx = ave_sx(1., Γ, β)
        writedlm(io,[β sx])
    end
end
