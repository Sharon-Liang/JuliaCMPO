using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using Printf
using PyPlot

z = make_operator(pauli(:z), 8)
β = 10.0
key = string(β)

g = 2
path = @sprintf "./data/g_%.1f.jld" g
d1 = load(path);
w = TFIsing(1.,g)
ψ = d1[key][2] |> tocmps

omega = [i for i in range(1.e-5,1,length=100)]

a = [spectral_function(i,z,z,ψ,w,β) for i in omega]
figure()
plot(omega, a, label = "spectrual function")
ta = @sprintf "spectral function : g = %.1f" g
title(ta)
legend()
PyPlot.display_figs()


s1 = [structure_factor(i,z,z,ψ,w,β,method =:S) for i in omega]
s2 = [structure_factor(i,z,z,ψ,w,β,method =:F) for i in omega]

figure()
plot(omega, s1, label = "S")
plot(omega, s2, "--", label = "F")
t = @sprintf "g = %.1f" g
title(t)
legend()
PyPlot.display_figs()


"""
    Nfreq = 20
    path1 = @sprintf "./data/Acmpo.txt"
    path2 = @sprintf "./data/input.txt" # frequency real_part imag_part \n, increasing masubara frenquancies

    open(path1, "w") do file
    end

    open(path2, "w") do file
    end

    for i in range(-10,10,step=0.05)
        open(path1, "a") do file1
            writedlm(file1,[i spectral_function(i,z,z,ψ,w,β,η=0.001)])
        end
    end

    for i = 1:Nfreq
        freq = Masubara_freq(i,β) |> imag
        gn = Masubara_GF(i,z,z,ψ,w,β)
        open(path2, "a") do file2
            writedlm(file2,[freq real(gn) imag(gn)])
        end
    end
"""