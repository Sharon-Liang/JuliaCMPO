using DelimitedFiles
using JLD, HDF5
using Printf


# rename file
model = "ising"
beta =[1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0]
gamma=[1.0, 1.5, 2.0]


folder = "imagtime"
func = "dReG"
pdir = @sprintf "./data/%s/%s/%s" model folder func
dir = @sprintf "./data/%s/%s/%s" model folder func
isdir(dir) || mkdir(dir)

D = 8
#eta = [0.001, 0.050]

op = [:+, :-, :x, :iy, :z]
nop= ["pm", "mp", "px", "py", "pz"]

for g = 1:length(gamma)
    for b = 1:length(beta), o = 1:length(op)
        β = beta[b]
        path1 = @sprintf "%s/g_%.1f_%s_D_%i_beta_%i.txt" pdir gamma[g] nop[o] D β
        path2 = @sprintf "%s/g_%.1f_%s_D_%im2_beta_%i.txt" pdir gamma[g] nop[o] D β 
        if op[o] == :z 
            npath1 = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" pdir gamma[g] D β 
            mv(path1, npath1)
            npath2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i.txt" pdir gamma[g] D β 
            mv(path2, npath2)
        else
            rm(path1)
            rm(path2)
        end
    end
end


for g = 1:length(gamma)
    (d0, r0) = divrem(g, 1)
    for d = 1:length(D)
        for b = 1:length(beta)
            β = beta[b]
            path1 = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir gamma[g] D[d] β
            data = readdlm(path1)
            data[:,1] = data[:,1] .* β
            open(path1, "w") do file
                for i = 1:length(data[:,1]) 
                    writedlm(file, [data[i,1] data[i,2]])
                end
            end
        end
    end
end












#delete files
#rm(path)