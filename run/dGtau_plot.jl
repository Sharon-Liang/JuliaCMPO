using DelimitedFiles
using JLD, HDF5
using Printf
using Plots 

D = 8

model = "ising"
folder = "imagtime"
dir = @sprintf "./data/%s/%s" model folder

g = 1.0
β = 10

path1 = @sprintf "%s/gtau/g_%.1f_D_%i_beta_%i.txt" dir g D β
gtau = readdlm(path1)

path2 = @sprintf "%s/dGtau/g_%.1f_D_%i_beta_%i.txt" dir g D β
dg = readdlm(path2)


len = 1601
dg1 = zeros(1600, 2)
for i = 1:len-1
    dg1[i,1] = gtau[i, 1]
    dtau = gtau[i+1, 1] - gtau[i, 1]
    dg1[i,2] = (gtau[i+1, 2] - gtau[i, 2])/dtau
end

plot(dg[:,1], dg[:,2], lw = 2, label=false)
plot!(dg1[:,1], dg1[:,2], lw = 2, label="finit diff")


plot(dg1[:,1], abs.(dg1[:,2] .- dg[1:len-1,2]), lw=2)
plot!(yaxis=:log)