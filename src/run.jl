using cMPO
using Optim
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using PyPlot

χ = 10

β = 20
#η = 1.e-4

len = 41
g = [i for i in range(0,2,length = len)]
#W = TFIsing(1.,1.)

#len = 100
#β = [i for i in range(0.1,10,length = len)]

#len = 50
#β = [i for i in range(0.01,20,length = len)]


open("./data/beta-20-Gamma.txt","w") do io
    jldopen("./data/beta-20-states.jld","w") do file
    ψ = init_cmps(χ)
    arr = toarray(ψ)
    for j in g
        W = TFIsing(1.0, j)
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

d = load("./data/beta-20-states.jld")
τ = [i for i in range(0, β, length = 100)]
key = string(0.10)
ψ = cmps(d[key][:,:,1],d[key][:,:,2])
W = TFIsing(1.0, 0.1)
open("./data/SzSz-beta-20-gamma-01.txt","w") do io
    for t in τ
        sz = Correlation_2time(pauli('z'), pauli('z'), ψ, W, β, t)
        writedlm(io,[t sz])
    end
end
