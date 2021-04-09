using cMPO
using DelimitedFiles
using HDF5, JLD
using PyPlot
using Printf

z = pauli('z')
β = 20.0 ; key = string(β)

W = TFIsing(1.0,1.0)
#load g=1.0 data
d1 = load("./data/g_1.0.jld")
ψ = cmps(d1[key][2][:,:,1],d1[key][2][:,:,2])

# theoretical C(0,t)
t0 = [x for x in range(0.2, β-0.2, step = 0.2)]
c0 = [critical_zz_cor(β,x,0) for x in t0]

# numerical C(0,t)
t1 = [x for x in range(0, β, step = 0.2)]
c1 = [correlation_2time(z,z,x,ψ,W,β) for x in t1]

figure()
plot(t0, c0, linewidth = 2, label = "theoretical")
plot(t1, c1, "--", linewidth = 2, label = "numerical")
xlim([0,β])
xlabel("τ")
ylabel("< z(0,τ) z(0,0) >")
title("β=20, Γ = J = 1.0, date : 2020-04-09")
legend()
PyPlot.display_figs()
