using cMPO
using PyPlot
using DelimitedFiles
using JLD, HDF5
using Printf
using Random; Random.seed!()
using LinearAlgebra

g = 2.
path = @sprintf "./data/g_%.1f.jld" g
d = load(path)
w = TFIsing(1.,g)
beta = [i for i in range(1,20,step=0.2)]
z = make_operator(pauli(:z),8)

path2 = @sprintf "./data/struc_fac_g_%.1f.txt" g
open(path2, "w") do file
end


R = 50 #random sampling times
D = 4 # matrix dimension
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
        open(path2, "a") do file
            writedlm(file,[β r S[1] S[2] S[3] S[4]])
        end
    end
end


d2 = readdlm(path2)
s0 = [π*b/4*correlation_2time(b/2,z,z,ψ,w,b) for b in beta]
figure()
scatter(1 ./ d2[:,1], d2[:,4])
#plot(1 ./ beta, s0, "-r",linewidth=2, label="β/2")
legend()
xlabel("T")
ylabel("S'(0)")
legend()
title(@sprintf "g = %.1f" g)
fname = @sprintf "struc_fac_g_%.1f_d1.pdf" g
PyPlot.savefig(fname,bbox_inches="tight")
PyPlot.display_figs()