include("./StructureFactorTrain.jl")
using DelimitedFiles, Printf
using NonNegLeastSquares
using LinearAlgebra

g = 1.0
β = 10.0
D = 8
ω = build_range(0.01, 10)

#load χ(τ) and dχ(τ) data
p1 = @sprintf "./data/ising/imagtime/gtau/g_%.1f_D_%i_beta_%i.txt" g D β
p2 = @sprintf "./data/ising/imagtime/dgtau/g_%.1f_D_%i_beta_%i.txt" g D β

d1 = readdlm(p1)
d2 = readdlm(p2)

x1 = d1[:,1]; y1 = d1[:,2]; data1 = [(x1, y1)]
x2 = d2[:,1]; y2 = d2[:,2]; data2 = [(x1, y1, x2, y2)]


#build kernel
dω = 0.01
ωmax = 100
ω = build_range(dω, ωmax)
len = length(ω)
K0 = build_kernal(0, x1, ω, β)
@time S = nonneg_lsq(K0,y1;alg=:nnls) 
y = K0 * S

#select α
eye = Matrix(1.0I, len, len)
ȳ1 = vcat(y1, zeros(len))

K̄ = vcat(K0, α*eye)
S̄ = nonneg_lsq(K̄,ȳ1;alg=:nnls) 
x2 = S̄' * S
ȳ = K0 * S̄
y2 = ȳ' * ȳ


