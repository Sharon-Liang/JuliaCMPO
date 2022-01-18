using DelimitedFiles, Printf
using Plots; gr(xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)

g = 1.0
β = 10.0
D = 8

pa = @sprintf "./spectral/data/Aw_g_%.1f_beta_%i.txt" g β
da = readdlm(pa)

ps = @sprintf "./spectral/data/Sw_g_%.1f_beta_%i.txt" g β
ds = readdlm(ps)

#read structure factor data
#chi optim
spc1 = @sprintf "./spectral/data/chi/g_%.1f_D_%i_beta_%i.txt" g D β
dc1 = readdlm(spc1)
Ac1 = ( 1 .- exp.(-β .* dc1[:,1])) .* dc1[:,2]

spc2 = @sprintf "./spectral/data/try1/g_%.1f_D_%i_beta_%i_chi.txt" g D β
dc2 = readdlm(spc2)
Ac2 = ( 1 .- exp.(-β .* dc2[:,1])) .* dc2[:,2]

#dchi optim
sp1 = @sprintf "./spectral/data/dchi/g_%.1f_D_%i_beta_%i.txt" g D β
d1 = readdlm(sp1)
A1 = ( 1 .- exp.(-β .* d1[:,1])) .* d1[:,2]

sp2 = @sprintf "./spectral/data/try1/g_%.1f_D_%i_beta_%i.txt" g D β
d2 = readdlm(sp2)
A2 = ( 1 .- exp.(-β .* d2[:,1])) .* d2[:,2]

plot(ds[:,1], ds[:,2], line=(:black, 2), label="Li")
plot!(dc1[:,1], dc1[:,2], line=(:dash, 2), label="χ")
plot!(dc2[:,1], dc2[:,2], line=(:dash, 2), label="χ")
plot!(d1[:,1], d1[:,2], line=(:solid, 1.5), label="χ+dχ")
plot!(d2[:,1], d2[:,2], line=(:solid, 1.5), label="χ+dχ")
plot!(xlabel="ω", ylabel="S(ω)", xlim=(0, 6))
fpath = @sprintf "./spectral/figure/Sw_g_%.1f_beta_%i.pdf" g β
Plots.savefig(fpath)

plot(da[:,1], da[:,2], line=(:black, 2), label="Li")
plot!(dc1[:,1], Ac1, line=(:dash, 2), label="χ")
plot!(dc2[:,1], Ac2, line=(:dash, 2), label="χ")
plot!(d1[:,1], A1, line=(:solid, 1.5), label="χ+dχ")
plot!(d2[:,1], A2, line=(:solid, 1.5), label="χ+dχ")
plot!(xlabel="ω", ylabel="A(ω)",xlim=(-0.1, 6))
fpath = @sprintf "./spectral/figure/Aw_g_%.1f_beta_%i.pdf" g β
Plots.savefig(fpath)