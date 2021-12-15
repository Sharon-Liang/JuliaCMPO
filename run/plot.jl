using Plots; pyplot(xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
default(palette = palette(:okabe_ito))
using DelimitedFiles
using HDF5, JLD
using Printf
using cMPO

invT = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0]
lt = length(invT)
D = 8
s0 = zeros(lt,7) #beta, li, Liang 8, 2m8
for t = 1:lt
    β = invT[t]
    p0 = @sprintf "./data/ising/spectrum-Li/Sw/g_1.0_beta_%i.txt" β
    p1 = @sprintf "./data/ising/spectrum-Li/Sw/g_1.0_beta_%i_eta_2pidbeta.txt" β
    p2 = @sprintf "./data/ising/imagtime/dReG/g_1.0_D_8_beta_%i.txt" β
    #p3 = @sprintf "./data/ising/imagtime/dReG_new/g_1.0_pz_D_8_beta_%i.txt" β

    #p4 = @sprintf "./data/ising/imagtime/gdivwn/g_1.0_D_8_beta_%i.txt" β
    #p5 = @sprintf "./data/ising/imagtime/gdivwn/g_1.0_D_8m2_beta_%i.txt" β

    s0[t,1] = 2π/β
    d0 = readdlm(p0); s0[t,2] = d0[800,2]
    d1 = readdlm(p1); s0[t,3] = d1[800,2]
    d2 = readdlm(p2); s0[t,4] = 2/β*d2[1,2]
    #d3 = readdlm(p3); s0[t,5] = 2/β*d3[1,2]
    #d4 = readdlm(p4); s0[t,6] = -2/β*d4[1,3]
    #d5 = readdlm(p5); s0[t,7] = -2/β*d5[1,3]
end

#scatter(s0[:,1], s0[:,2], label="η=0.001")
scatter(s0[:,1], s0[:,3], marker=(6,stroke(0.5)),label="η=2π/β")
scatter!(s0[:,1], s0[:,4], marker=(6,stroke(0.5)), label="dReG D=8")
scatter!(s0[:,1], s0[:,5], marker=(:star, 6,stroke(0.5)), label="dReG_new D=8")

#scatter!(s0[:,1], s0[:,6], marker=(6,stroke(0.5)), label="G/ΔE D=8")
#scatter!(s0[:,1], s0[:,7], marker=(:star, 6,stroke(0.5)), label="G/ΔE D=8×2")

plot!(xlabel="2π/β", ylabel="S(0)")
savefig("./figure/compare_S0_TdRe_wrong_dRe.pdf")


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


g = 1.0
w = TFIsing(1.0, g)
data_path = @sprintf "./data/ising/D_%i/g_%.1f.jld" D g
data = load(data_path)
st = similar(invT)
z = make_operator(pauli(:z),D);
for t = 1:lt
    β = invT[t]
    η = 0.001
    β = invT[t]; key = string(β)
    ψ1 = tocmps(data[key][2])
    st[t] = structure_factor(0,z',z, ψ1, w, β, η = η)
end

scatter(s0[:,1], s0[:,2],marker=(6,stroke(0.5)),label="Li, η=0.001")
scatter!(s0[:,1], st,marker=(:star, 6,stroke(0.5)),label="cmpo, η=0.001")
plot!(xlabel="2π/β", ylabel="S(0)")
savefig("./figure/compare_S0_TdRe_S0_eta_0p001.pdf")


delta(0, 0.001)
delta(0,2π)