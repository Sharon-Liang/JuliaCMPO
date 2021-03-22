using LinearAlgebra
using StatsFuns
using DelimitedFiles
using JLD, HDF5

include("exact.jl")

# exact free energy
β = 20
len = 51
g = [i for i in range(0,2,length = len)]

open("./data/0323/fexact-beta-20.txt","w") do io
    for j in g
        f = free_energy(1., j, β)
        writedlm(io,[j f])
    end
end


# exact sx
β = 20
len = 51
g = [i for i in range(0,2,length = len)]

open("./data/0323/sx-exact-beta-20.txt","w") do io
    for j in g
        sx = ave_sx(1., j, β)
        writedlm(io,[j sx])
    end
end
