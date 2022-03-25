using LinearAlgebra
using DelimitedFiles, Printf, Plots
include("/home/sliang/JuliaCode/mycMPO/src/KernelFunction.jl")

g = 1.0
β = 10.0
const D = 8
const dω = 0.01
const ωmax = 10.0
const DataFolder = "/data/sliang/CMPO/ising/imagtime"
const FigureFolder = "/home/sliang/JuliaCode/mycMPO/figure_nnls"
#DataFile Path 
DataFile1 = @sprintf "%s/gtau/g_%.1f_D_%i_beta_%i.txt" DataFolder g D β

#LOAD DATA
d1 = readdlm(DataFile1); x1 = d1[:,1]; y1 = d1[:,2]

#BUILD KERNEL
#ω, len = build_range(dω, ωmax)
#K1 = build_kernal(0, x1, ω, β)
ω, len = build_range(dω, ωmax, sqrt(eps()))
K1 = build_Akernal(x1, ω, β)

#Piscard Plot
F1 = svd(K1)  #Store the Factorization Object
ub = abs.(transpose(F1.U) * y1)
ub_s = ub ./ abs.(F1.S) 
id = [i for i=1:length(ω)]
figtitle = @sprintf "Piscard plot: g=%.1f, D=%i, β=%i" g D β
plot(id, F1.S, line=(1)
    ,marker=(:dot)
    ,label="σi")

plot!(id, ub, line=(1)
    ,marker=(:cross)
    ,label="|ui^T b|")

plot!(id, ub_s, line=(1)
    ,marker=(:star)
    ,label="|ui^T b/σi|")

plot!(xlabel="i", xlim=(1,30)
    ,yaxis=:log
    ,title=figtitle
    ,framestyle=:box
    ,minorgrids=true
    ,legend=:bottomleft)

FigureName = @sprintf "ising_piscard_g_%.1f_D_%i_beta_%i_Akernel.pdf" g D β
savefig(FigureFolder*"/"*FigureName)