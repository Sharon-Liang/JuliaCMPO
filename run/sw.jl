include("../src/StructureFactorTrain.jl")

println("2022.01.11: TFIsing: structure factor train")
g = 1.0
β = 10.0
D = 8

#load χ(τ) and dχ(τ) data
p1 = @sprintf "../data/ising/imagtime/gtau/g_%.1f_D_%i_beta_%i.txt" g D β
p2 = @sprintf "../data/ising/imagtime/dgtau/g_%.1f_D_%i_beta_%i.txt" g D β

d1 = readdlm(p1)
d2 = readdlm(p2)

x1 = d1[:,1]; y1 = d1[:,2]; data1 = [(x1, y1)]
x2 = d2[:,1]; y2 = d2[:,2]; data2 = [(x1, y1, x2, y2)]

dir = "../data/ising/spectrum_nn"
isdir(dir) || mkdir(dir)

model, parameters = build_NN_model()
s(ω::Real) = model([ω])[1]

err = loss(data1[1],β,s)
ep = @sprintf "%s/g_%.1f_D_%i_beta_%i_error.txt" dir g D β
jldopen(ep, "w") do file
    writedlm(file, err)
end

tolerance = 1.e-5
step = 100
loop = 100
epoch = 0
for l = 1:loop, i = 1:step
    epoch += 1
    train!((x,y) -> loss((x,y),β,s), parameters, data1, Descent(0.1))
    err = loss(data1[1],β,s)
    println(epoch, ": ", err )
    jldopen(ep, "w+") do file
        writedlm(file, err)
    end
    if err < tolerance
        break
    end
end


