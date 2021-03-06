using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5

β = 20
len = 51
g = [i for i in range(0,2,length = len)]

# Thermal_average
d = load("./data/beta-20-chi-8/srandom-z.jld")
open("./data/beta-20-chi-8/sz-rand.txt","w") do io
    for j in g
        key = string(j)
        ψ = cmps(d[key][:,:,1],d[key][:,:,2])
        W = TFIsing(1.0, j,field ='z')
        sz = Thermal_average(ψ,W,pauli('z'),β)
        writedlm(io,[j sz])
    end
end

# 2-time correlations
d = load("./data/beta-20-chi-8/srandom.jld")
tau = [i for i in range(0, β, length = 100)]
open("./data/beta-20-chi-8/zz-rand-2.txt","w") do io
    for τ in tau
        j = g[50]
        key = string(j)
        ψ = cmps(d[key][:,:,1],d[key][:,:,2])
        W = TFIsing(1.0, j)
        corr = Correlation_2time(pauli('z'), pauli('z'), ψ, W, β, τ)
        writedlm(io,[τ corr])
    end
end
