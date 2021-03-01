using cMPO
using Optim
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using PyPlot

χ = 10

β = 1.e4
η = 1.e-4

len = 41
g = [i for i in range(0,2,length = len)]
#W = TFIsing(1.,1.)

#len = 100
#β = [i for i in range(0.1,10,length = len)]

#len = 50
#β = [i for i in range(0.01,20,length = len)]


open("./data/zero-T-Gamma.txt","w") do io
    jldopen("./data/zero-T-states.jld","w") do file
    ψ = init_cmps(χ)
    arr = toarray(ψ)
    for j in g
        W = TFIsing(1.0, j, field = η)
        of(x::Array{Float64, 3}) = OptimFreeEnergy(x::Array{Float64, 3}, W, β)
        of!(gx::Array{Float64, 3}, x::Array{Float64,3}) = OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, β)
        op = optimize(of, of!, arr, LBFGS())
        writedlm(io,[j minimum(op)])
        arr = op.minimizer
        key = string(j)
        write(file, key, arr)
    end
    end
end


d = load("./data/zero-T-states.jld")
open("./data/Sz-Sz-tau.txt","w") do io
    for j in g
        key = string(g[j])
        ψ = cmps(d[key][:,:,1],d[key][:,:,2])
        W = TFIsing(1.0, g[j], field = η)
        sz = Thermal_average(ψ, W, pauli('z'), β)
        writedlm(io,[j sz])
    end
end
