using cMPO
using DelimitedFiles
using HDF5, JLD
using PyPlot
using Printf

z = pauli('z')

#load g=1.0, 0.1, 2.0 data
d1 = load("./data/g_1.0.jld")
d2 = load("./data/g_0.1.jld")
d3 = load("./data/g_2.0.jld")

J = 1.0; Γ = 1.0; W = TFIsing(J,Γ)

β = 20.0 ; key = string(β)
ψ = cmps(d1[key][2][:,:,1], d1[key][2][:,:,2])

w = [x for x in range(1.e-3,10,length = 1000)]
c = [critical_zz_chi(x,β) for x in w]
# theoretical C(0,t)
c0 = [imag_susceptibility(x,z,z,ψ,W,β,η = 0.05) for x in w]
# numerical C(0,t)
c1 = [imag_susceptibility(x,z,z,ψ,W,β,η = 0.1) for x in w]

figure()
semilogx(w, c, "k", label = "theoretical")
semilogx(w, c0, "o",markersize = 3, label = "η = 0.05")
semilogx(w, c1, "o",markersize = 3, label = "η = 0.1")
xlabel("ω")
ylabel("Im χ(ω)")
title("β = 20, Γ = J = 1.0 , date : 2020-04-12")
legend()
PyPlot.display_figs()

figure()

xlabel("ω")
ylabel("Im χ(ω)")
title("β = 20, Γ = J = 1.0 , date : 2020-04-12")
legend()
PyPlot.display_figs()
