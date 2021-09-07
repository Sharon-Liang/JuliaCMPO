using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using Printf
#using PyPlot

g = 2.0; w = TFIsing(1.,g)
β = 10.0; key = string(β)

z8 = make_operator(pauli(:z), 8)
z16 = make_operator(pauli(:z), 16)

path1 = @sprintf "./data/g_%.1f.jld" g
path2 = @sprintf "./data/chi16/g_%.1f.jld" g

d8 = load(path1)
d16 = load(path2)

ψ8 = d8[key][2] |> tocmps
ψ28 = w * ψ8
ψ16 = d16[key][2] |> tocmps

#structure_factor
omega = [10^i for i in range(-5,1,length=100)]
s8 = [structure_factor(i,z8,z8,ψ8,w,β,method =:S) for i in omega]
s28 = [structure_factor(i,z16,z16,ψ28,w,β,method =:S) for i in omega]
s16 = [structure_factor(i,z16,z16,ψ16,w,β,method =:S) for i in omega]

path3 = @sprintf "./data/Sw_g_%.1f.txt" g
open(path3, "w") do file
    for i = 1:length(omega)
        writedlm(file,[omega[i] s8[i] s28[i] s16[i]])
    end
end

# Correlation_2time
omega = [i for i in range(0,β,length=100)]
s8 = [correlation_2time(i,z8,z8,ψ8,w,β) for i in omega]
s28 = [correlation_2time(i,z16,z16,ψ28,w,β) for i in omega]
s16 = [correlation_2time(i,z16,z16,ψ16,w,β) for i in omega]

path3 = @sprintf "./data/Gt_g_%.1f.txt" g
open(path3, "w") do file
    for i = 1:length(omega)
        writedlm(file,[omega[i] s8[i] s28[i] s16[i]])
    end
end
