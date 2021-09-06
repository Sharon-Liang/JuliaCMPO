using cMPO
using PyPlot
using DelimitedFiles
using JLD, HDF5
using Printf

z = make_operator(pauli(:z),8)
β = 10
path = @sprintf "./data/b_%i_jchange.jld" β
d = load(path)

J = 0.; key =string(J)

ψ = tocmps(d[key][2])
w = TFIsing(0.,1.)

tau = [t for t in range(0,β,length=100)]

gt1 = [correlation_2time(t,z,z,ψ,w,β) for t in tau]
gt2 = [correlation_2time(0,t,β) for t in tau]

dg = gt1 .- gt2

ψ0 = init_cmps(2,w)
z0 = make_operator(pauli(:z),ψ0)
gt3 = [correlation_2time(t,z0,z0,ψ0,w,β) for t in tau]

dg2 = gt3 .- gt2

τ = tau ./ β
figure()
plot(τ, gt2, linewidth=2, label="theory")
plot(τ, gt1, "--o", linewidth = 2, label = "χ=8 cmpo")
plot(τ, gt3, "--*", linewidth = 2, label = "χ=2 cmpo")
xlabel("τ/β")
ylabel("G(τ)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_gf.pdf",bbox_inches="tight")
PyPlot.display_figs()

figure()
plot(τ, abs.(dg), linewidth = 2, label = "χ=8 cmpo")
plot(τ, abs.(dg2),linewidth = 2, label = "χ=2 cmpo")
xlabel("τ/β")
ylabel("G(τ) difference (abs)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_gf_diff.pdf",bbox_inches="tight")
PyPlot.display_figs()

#spectral_function

omega = [i for i in range(-4,4, length=100)]
a1 = [spectral_function(0.,x,β) for x in omega]
a2 = [spectral_function(x,z,z, ψ, w,β) for x in omega]
a3 = [spectral_function(x,z0,z0, ψ0, w,β) for x in omega]

figure()
plot(omega, a1, linewidth=2, label="theory")
plot(omega, a2, "--o", linewidth = 2, label = "χ=8 cmpo")
plot(omega, a3, "--*", linewidth = 2, label = "χ=2 cmpo")
xlabel("ω")
ylabel("Imχ(ω)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_imag_chi.pdf",bbox_inches="tight")

figure()
plot(omega, abs.(a2 .- a1), linewidth = 2, label = "χ=8 cmpo")
plot(omega, abs.(a3 .- a1),"--", linewidth = 2, label = "χ=2 cmpo")
xlabel("ω")
ylabel("Imχ(ω) difference (abs)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_imag_chi_diff.pdf",bbox_inches="tight")

