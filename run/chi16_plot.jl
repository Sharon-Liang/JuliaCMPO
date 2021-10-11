using cMPO
using PyPlot
using DelimitedFiles
using JLD, HDF5
using Printf

β = 20
g = 1.0
type = "Gt"
if type == "Sw"
    xlab = "ω" ; ylab = "S(ω)"; 
elseif type == "Gt"
    xlab = "τ/β" ; ylab = "χ(τ)"
elseif type == "floss"
    xlab = "β" ; ylab = "floss"
elseif type == "ImX"
    xlab = "ω" ; ylab = "Imχ(ω)"
end

path = @sprintf "./data/%s_g_%.1f_β_%i.txt" type g β
d = readdlm(path)

figure()
plot(d[:,1]./β, d[:,2], linewidth=3,label = "χ=8")
##plot(d[:,1], d[:,3], linewidth=2, label = "χ=2×8")
#plot(d[:,1], d[:,4], linewidth=2, label = "χ=16")
#exact results
if g == 1. && type == "ImX"
    omega = [i for i in range(-5,5,length=200)]
    s0 = [spectral_density(g,i,β) for i in omega]
    plot(omega, s0, "--k", linewidth=2, label="theory")
elseif g == 1. && type == "Gt"
    tau = [i for i in range(0.1,β-0.1,length=100)]
    g0 = [correlation_2time(g,i,β) for i in tau]
    plot(tau ./β, g0, "--k", linewidth=3, label="theory")
end

xlabel(xlab)
ylabel(ylab)
#xlim(0,5)
#ylim(-0.1,)
type == "floss" ? (tit=@sprintf "g = %.1f" g) : (tit = @sprintf "g = %.1f, β = %i" g β)
#type == "Sw" ? xscale("log") : 
title(tit)
legend()
fname = @sprintf "%s_g_%.1f_pos.pdf" type g
PyPlot.savefig(fname,bbox_inches="tight")
PyPlot.display_figs()