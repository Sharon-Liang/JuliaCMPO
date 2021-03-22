using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5

Γ = 2.0
len = 40
beta = [i for i in range(1,20,length = len)]

# Thermal_average
d = load("./data/0323/gamma-2.0.jld")
open("./data/0323/sx-gamma-2.0.txt","w") do io
    for β in beta
        key = string(β)
        ψ = cmps(d[key][:,:,1],d[key][:,:,2])
        W = TFIsing(1.0, Γ)
        sx = Thermal_average(ψ,W,pauli('x'),β)
        writedlm(io,[β sx])
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
