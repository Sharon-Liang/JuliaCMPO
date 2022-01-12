using DelimitedFiles, Printf
include("../src/StructureFactorTrain.jl")


println("2022.01.11: TFIsing: structure factor train")
g = 1.0
β = 10.0
D = 8
ω = build_range(0.01, 10)

#load χ(τ) and dχ(τ) data
p1 = @sprintf "../data/ising/imagtime/gtau/g_%.1f_D_%i_beta_%i.txt" g D β
p2 = @sprintf "../data/ising/imagtime/dgtau/g_%.1f_D_%i_beta_%i.txt" g D β

d1 = readdlm(p1)
d2 = readdlm(p2)

x1 = d1[:,1]; y1 = d1[:,2]; data1 = [(x1, y1)]
x2 = d2[:,1]; y2 = d2[:,2]; data2 = [(x1, y1, x2, y2)]

dir = "../data/ising/spectrum_nn"
isdir(dir) || mkdir(dir)

dir1 = @sprintf "%s/chi1" dir
isdir(dir1) || mkdir(dir1)

c = 0
err = 100
while err > 1.0
    c += 1
    model, parameterr = build_NN_model()
    s(ω::Real) = model([ω])[1]
    err = loss(data2[1],β,s)
end
println("try ", cont, " time to initiate the model")

sp = @sprintf "%s/g_%.1f_D_%i_beta_%i_rand.txt" dir1 g D β
open(sp, "w") do file
    for i = 1:length(ω)
        writedlm(file, s(ω[i]))
    end
end


ep = @sprintf "%s/g_%.1f_D_%i_beta_%i_error.txt" dir1 g D β
open(ep, "w") do file
    writedlm(file, err)
end

tolerance = 1.e-5
epoch_max = 10000
for epoch = 1: epoch_max
    train!((x1,y1, x2, y2) -> loss((x1,y1, x2, y2),β,s), parameterr, data1, Descent())
    err = loss(data2[1],β,s)
    println(epoch, ": ", err )
    open(ep, "a") do file
        writedlm(file, err)
    end
    if err < tolerance
        break
    end
end


sp = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir1 g D β
open(sp, "w") do file
    for i = 1:length(ω)
        writedlm(file, s(ω[i]))
    end
end


using BSON: @save
mp = @sprintf "%s/g_%.1f_D_%i_beta_%i.bson" dir1 g D β
@save mp model
