using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using Printf
#using PyPlot
z8 = make_operator(pauli(:z), 8)
z16 = make_operator(pauli(:z), 16)
gamma = [0.5, 1.0, 2.0]
for g in gamma
    w = TFIsing(1.,g)
    β = 10.0; key = string(β)

    path1 = @sprintf "./data/g_%.1f.jld" g
    path2 = @sprintf "./data/chi16/g_%.1f.jld" g
    d8 = load(path1)
    d16 = load(path2)

    ψ8 = d8[key][2] |> tocmps
    ψ28 = w * ψ8
    ψ16 = d16[key][2] |> tocmps

    #structure_factor
    omega = [i for i in range(0,5,length=100)]
    s8 = [structure_factor(i,z8,z8,ψ8,w,β,method =:S) for i in omega]
    s28 = [structure_factor(i,z16,z16,ψ28,w,β,method =:S) for i in omega]
    s16 = [structure_factor(i,z16,z16,ψ16,w,β,method =:S) for i in omega]

    path3 = @sprintf "./data/Sw_g_%.1f.txt" g
    open(path3, "w") do file
        for i = 1:length(omega)
            writedlm(file,[omega[i] s8[i] s28[i] s16[i]])
        end
    end

    #spectral density
    omega = [i for i in range(-5,5,length=200)]
    s8 = [spectral_density(i,z8,z8,ψ8,w,β) for i in omega]
    s28 = [spectral_density(i,z16,z16,ψ28,w,β) for i in omega]
    s16 = [spectral_density(i,z16,z16,ψ16,w,β) for i in omega]

    path3 = @sprintf "./data/ImX_g_%.1f.txt" g
    open(path3, "w") do file
        for i = 1:length(omega)
            writedlm(file,[omega[i] s8[i] s28[i] s16[i]])
        end
    end
end

"""
# free_energy
omega = [i for i in range(1,20,length=100)]
fe = [free_energy(1.,g,i) for i in omega]
s8 = [free_energy(ψ8,w,i) for i in omega] .- fe
s28 = [free_energy(ψ28,w,i) for i in omega] .- fe
s16 = [free_energy(ψ16,w,i) for i in omega] .- fe

path3 = @sprintf "./data/floss_g_%.1f.txt" g
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
"""

