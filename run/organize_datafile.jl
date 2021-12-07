using DelimitedFiles
using JLD, HDF5
using Printf




# rename file

model = "ising"
D = [8]
beta =[10, 20, 30, 40]
gamma=[1.0, 1.5, 2.0]

folder = "imagtime-hessian"
func = "gtau"
#pdir = @sprintf "./data/%s/%s" model folder
dir = @sprintf "./data/%s/%s/%s" model folder func
#isdir(dir) || mkdir(dir)


for g = 1:length(gamma)
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









for g = 1:length(gamma)
    for d = 1:length(D)
        for b = 1:length(beta)
            β = beta[b]
            path1 = @sprintf "%s/%s_g_%.1f_D_%i_beta_%i.txt" pdir func gamma[g] D[d] β
            npath1 = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir gamma[g] D[d] β
            mv(path1, npath1)
            path2 = @sprintf "%s/%s_g_%.1f_D_2m%i_beta_%i.txt" pdir func gamma[g] D[d] β
            npath2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i.txt" dir gamma[g] D[d] β
            mv(path2, npath2)
        end
    end
end


#delete files
#rm(path)