using cMPO
using PyPlot
using DelimitedFiles
using JLD, HDF5
using Printf

β = 10
g = 0.5
type = "ImX"
if type == "Sw"
    xlab = "ω" ; ylab = "S(ω)"; 
elseif type == "Gt"
    xlab = "τ" ; ylab = "G(τ)"
elseif type == "floss"
    xlab = "β" ; ylab = "floss"
elseif type == "ImX"
    xlab = "ω" ; ylab = "Imχ(ω)"
end

path = @sprintf "./data/%s_g_%.1f.txt" type g
d = readdlm(path)

figure()
plot(d[:,1], d[:,2], linewidth=2, label = "χ=8")
plot(d[:,1], d[:,3], linewidth=2, label = "χ=2×8")
plot(d[:,1], d[:,4], linewidth=2, label = "χ=16")
#exact results
if g == 1. && type == "ImX"
    omega = [i for i in range(-5,5,length=200)]
    s0 = [spectral_density(g,i,β) for i in omega]
    plot(omega, s0, "--k", linewidth=2, label="theory")
elseif g == 1. && type == "Gt"
    tau = [i for i in range(0.1,β-0.1,length=100)]
    g0 = [correlation_2time(g,i,β) for i in tau]
    plot(tau, g0, "--k", linewidth=2, label="theory")
end

xlabel(xlab)
ylabel(ylab)
xlim(0,5)
ylim(-0.1,)
type == "floss" ? (tit=@sprintf "g = %.1f" g) : (tit = @sprintf "g = %.1f, β = 10" g)
#type == "Sw" ? xscale("log") : 
title(tit)
legend()
fname = @sprintf "%s_g_%.1f_pos.pdf" type g
PyPlot.savefig(fname,bbox_inches="tight")
PyPlot.display_figs()