using cMPO
using Optim
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using PyPlot

χ = 10

β = 20
#η = 1.e-4

len = 50
g = [i for i in range(0,2,length = len)]
#W = TFIsing(1.,1.)

#len = 100
#β = [i for i in range(0.2,20,length = len)]

#len = 50
#β = [i for i in range(0.01,20,length = len)]


#open("./data/Gamma-1/F_energy.txt","w") do io
jldopen("./data/beta-20/states-xfield.jld","w") do file
    ψ = init_cmps(χ)
    arr = toarray(ψ)
    for j in g
        W = TFIsing(1.0, j, field = 'x')
        of(x::Array{Float64, 3}) = OptimFreeEnergy(x::Array{Float64, 3}, W, β)
        of!(gx::Array{Float64, 3}, x::Array{Float64,3}) = OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, β)
        op = optimize(of, of!, arr, LBFGS())
        #writedlm(io,[b minimum(op)])
        arr = op.minimizer
        key = string(j)
        write(file, key, arr)
    end
    #end
end


d = load("./data/beta-20/states-xfield.jld")
#τ = [i for i in range(0, β, length = 100)]
open("./data/beta-20/Sx.txt","w") do io
    for j in g
        key = string(j)
        ψ = cmps(d[key][:,:,1],d[key][:,:,2])
        W = TFIsing(1.0, j, field='x')
        sz = Thermal_average(ψ, W, pauli('x'),β)
        writedlm(io,[j sz])
    end
end
