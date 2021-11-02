using cMPO
using PyPlot
using DelimitedFiles
using JLD, HDF5
using Printf
using Random; Random.seed!()
using LinearAlgebra

g = 1.0
beta = [10.0, 20.0, 30.0, 40.0]
path = @sprintf "./data/ug_%.1f.jld" g
d = load(path)
w = TFIsing(1.,g)
#beta = [i for i in range(1,20,step=0.2)]
z = make_operator(pauli(:z),16)

for i = 1:4
    β = beta[i]; key=string(β)
    ψ = d[key][2] |>tocmps
    ψ = w * ψ
    c = -Masubara_freq_GF(0, z,z,ψ, w, β)
    s = @sprintf "β = %i, c = %.5f" β c*2π
    println(s)
end

R = 50 #random sampling times
D = 8 # matrix dimension

path2 = @sprintf "./data/struc_fac_g_%.1f_order_%i.txt" g D
open(path2, "w") do file
    for i = 1:length(beta)
        β = beta[i]; key = string(β)
        ψ = tocmps(d[key][2])
        for r = 1:R
            τ = rand(D)
            while prod(τ) == 0
                τ = rand(D)
            end
            Gt = [π*correlation_2time(β*t,z,z,ψ,w,β) for t in τ]
            M = zeros(D,D)
            for i = 1:D, j = 1:D
                M[i,j] = β^(-j) *(τ[i]^(-j)+ (1-τ[i])^(-j))
            end
            S = M^(-1) * Gt
            writedlm(file,[β r S[1] S[2] S[3] S[4]])
        end
    end
end

s0 = zeros(length(beta))
for i = 1:length(beta)
    β = beta[i]; key = string(β)
    ψ = tocmps(d[key][2])
    s0[i] = π*β/4*correlation_2time(β/2,z,z,ψ,w,β)
end

d2 = readdlm(path2)
figure()
scatter(1 ./ d2[:,1], d2[:,3])
plot(1 ./ beta, s0, "-r",linewidth=2, label="β/2")
xlabel("T")
ylabel("S(0)")
ylim(-0.5, 0.5)
legend()
title(@sprintf "g = %.1f, order = %i" g D)
fname = @sprintf "struc_fac_g_%.1f_order_%i.pdf" g D
PyPlot.savefig(fname,bbox_inches="tight")
PyPlot.display_figs()

figure()
scatter(1 ./ d2[:,1], d2[:,4])
#plot(1 ./ beta, s0, "-r",linewidth=2, label="β/2")
xlabel("T")
ylabel("S'(0)")
legend()
title(@sprintf "g = %.1f" g)
fname = @sprintf "struc_fac_g_%.1f_d1_order_%i.pdf" g D
PyPlot.savefig(fname,bbox_inches="tight")
PyPlot.display_figs()
