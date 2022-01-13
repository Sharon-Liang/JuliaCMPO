include("./StructureFactorTrain.jl")
using DelimitedFiles, Printf
using Optim, Zygote
using FluxOptTools

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

#build NN model
model, pars = build_NN_model(5)
s(ω::Real) = model([ω])[1]

Zygote.refresh() # currently needed when defining new adjoints
lossfun, gradfun, fg!, p0 = optfuns(()->loss(data1[1], β, s), pars)
res = Optim.optimize(Optim.only_fg!(fg!), p0, LBFGS(), Optim.Options(iterations=10000, store_trace=true))
