using DelimitedFiles, Printf
using BSON: @save
include("./StructureFactor.jl")

println("2022.01.13: TFIsing: structure factor train")
g = 1.0
β = 10.0
D = 8
ω, _ = build_range(0.01, 10)

#load χ(τ) and dχ(τ) data
p1 = @sprintf "./data/ising/imagtime/gtau/g_%.1f_D_%i_beta_%i.txt" g D β
p2 = @sprintf "./data/ising/imagtime/dgtau/g_%.1f_D_%i_beta_%i.txt" g D β

d1 = readdlm(p1)
d2 = readdlm(p2)

x1 = d1[:,1]; y1 = d1[:,2]; data1 = [(x1, y1)]
x2 = d2[:,1]; y2 = d2[:,2]; data2 = [(x1, y1, x2, y2)]

dir = "../data/ising/spectrum_nn/dchi"
isdir(dir) || mkdir(dir)

opt = ADAM()
err = []
model, pars = build_NN_model(5)
s(ω::Real) = model([ω])[1]

#save the init model
mpi = @sprintf "%s/g_%.1f_D_%i_beta_%i_init.bson" dir g D β 
@save mpi model

#save init spectrum
spi = @sprintf "%s/sw_g_%.1f_D_%i_beta_%i_init.txt" dir g D β
open(spi, "w") do file
    for i = 1:length(ω)
        writedlm(file, [ω[i]  s(ω[i])])
    end
end
    
#creat error file
err = loss(data2,β,s)[1]
ep = @sprintf "%s/err_g_%.1f_D_%i_beta_%i.txt" dir g D β
open(ep, "w") do file
    writedlm(file, err)
end

tolerance = 1.e-5
epoch_max = 2000
loop = 10
#for l = 1:loop, epoch = 1: epoch_max
    train!((x,y,x1,y1) -> loss((x,y,x1,y1),β,s), pars, data2, opt)
    err = loss(data2[1],β,s)
    #println(epoch, ": ", err )
    open(ep, "a") do file
        writedlm(file, err)
    end
    if err < tolerance
        break
    end
#end

sp = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir g D β
open(sp, "w") do file
    for i = 1:length(ω)
        writedlm(file, s(ω[i]))
    end
end

mp = @sprintf "%s/g_%.1f_D_%i_beta_%i.bson" dir g D β
@save mp model
