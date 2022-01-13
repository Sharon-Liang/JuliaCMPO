include("./StructureFactorTrain.jl")
using DelimitedFiles, Printf
using Optim, Zygote
using FluxOptTools
using BSON: @save, @load

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

loss(data1[1],β,s)

dir = "./spectral/data/chi"
isdir(dir) || mkdir(dir)
#save the init model
mpi = @sprintf "%s/LBFGS_g_%.1f_D_%i_beta_%i_init.bson" dir g D β 
@save mpi model

Zygote.refresh() # currently needed when defining new adjoints
@time lossfun, gradfun, fg!, p0 = optfuns(()->loss(data1[1], β, s), pars)
@time res = Optim.optimize(Optim.only_fg!(fg!), p0, 
            Optim.Options(iterations=3000, store_trace=true, show_trace=true))

mp = @sprintf "%s/LBFGS_g_%.1f_D_%i_beta_%i.bson" dir g D β 
@save mp model
sp = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir g D β
open(sp, "w") do file
    for i = 1:length(ω)
        writedlm(file, [ω[i] s(ω[i])])
    end
end

trace = Optim.f_trace(res)
ep = @sprintf "%s/ftrace_g_%.1f_D_%i_beta_%i.txt" dir g D β
open(ep, "w") do file
    writedlm(file, trace)
end

#build NN model: dchi
model, pars = build_NN_model(5)
s(ω::Real) = model([ω])[1]
loss(data2[1],β,s)

dir = "./spectral/data/dchi"
isdir(dir) || mkdir(dir)

#save the init model
mpi = @sprintf "%s/LBFGS_g_%.1f_D_%i_beta_%i_init.bson" dir g D β 
@save mpi model


Zygote.refresh() # currently needed when defining new adjoints
@time lossfun, gradfun, fg!, p0 = optfuns(()->loss(data2[1], β, s), pars)
@time res = Optim.optimize(Optim.only_fg!(fg!), p0, LBFGS(),
            Optim.Options(iterations=3000, store_trace=true, show_trace=true))

mp = @sprintf "%s/LBFGS_g_%.1f_D_%i_beta_%i.bson" dir g D β 
@save mp model
sp = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir g D β
open(sp, "w") do file
    for i = 1:length(ω)
        writedlm(file, [ω[i] s(ω[i])])
    end
end
            
trace = Optim.f_trace(res)
ep = @sprintf "%s/ftrace_g_%.1f_D_%i_beta_%i.txt" dir g D β
open(ep, "w") do file
    writedlm(file, trace)
end