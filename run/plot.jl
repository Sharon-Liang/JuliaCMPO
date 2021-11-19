using Plots; pyplot(xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
default(palette = palette(:okabe_ito))
using DelimitedFiles
using Printf
using cMPO

g = 1.0
path = @sprintf "./data/floss_g_%.1f_new.txt" g
d = readdlm(path)

plot(d[:,1], abs.(d[:,2]),lw=2, marker=(:circle,3,stroke(0.1)),label="χ=8")
plot!(d[:,1], abs.(d[:,3]),lw =2, label="χ=8*2")
plot!(d[:,1], abs.(d[:,4]),lw=2, marker=(:circle,3,stroke(0.1)),label="χ=16")
plot!(d[:,1], abs.(d[:,5]),lw=2, label="χ=16*2")
plot!(d[:,1], abs.(d[:,6]),lw=2, label="χ=8*2×2")
plot!(d[:,1], abs.(d[:,7]),line=(:dash,2), marker=(:circle,3,stroke(0.1)),label="hessian χ=8")
plot!(d[:,1], abs.(d[:,8]),line=(:dash,2), label="hessian χ=8*2")
plot!(yscale=:log10)
plot!(xlabel="β", ylabel="free energy loss",legend=:bottomright)
ptitle=@sprintf "g = %.1f" g
plot!(title=ptitle)
pname =@sprintf "./notes/floss_g_%.1f.pdf" g
savefig(pname)