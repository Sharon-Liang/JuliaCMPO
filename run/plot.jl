using Plots; pyplot(xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=12)
default(palette = palette(:okabe_ito))
using DelimitedFiles
using Printf
using cMPO

g = 2.0 
path = @sprintf "./data/floss_g_%.1f.txt" g
d = readdlm(path)

plot(d[:,1], d[:,2],lw=2, label="χ=8")
plot!(d[:,1], d[:,3],lw =2, label="χ=8*2")
plot!(d[:,1], d[:,4],lw=2, label="χ=16")
plot!(d[:,1], d[:,5],lw=2, label="χ=16*2")
plot!(d[:,1], d[:,6],lw=2, label="χ=8*2×2")
plot!(lengend=:bottomright)
plot!(xlabel="β", ylabel="free energy loss")
plot!(ylim=(1.e-16,1.e-5), yscale=:log10)
