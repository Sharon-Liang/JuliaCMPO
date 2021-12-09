using DelimitedFiles
using JLD, HDF5
using Printf


# rename file
model = "ising"
beta =[10, 20, 30, 40]
gamma=[1.0, 1.5, 2.0]


folder = "spectrum-direct-ac"
func = "Sw"
pdir = @sprintf "./data/%s/%s" model folder
dir = @sprintf "./data/%s/%s/%s" model folder func
isdir(dir) || mkdir(dir)

eta = [0.001, 0.050]
D = 8
for g = 1:length(gamma)
    (d0, r0) = divrem(gamma[g], 1)
    for b = 1:length(beta), e = 1:length(eta)
        β = beta[b]
        η = eta[e]
        path1 = @sprintf "%s/%s_g_%.1f_D_%i_beta_%i_eta_%.3f.txt" pdir func gamma[g] D β η
        npath1 = @sprintf "%s/g_%.1f_D_%i_beta_%i_eta_%.3f.txt" dir gamma[g] D β η
        mv(path1, npath1)
        path2 = @sprintf "%s/%s_g_%.1f_D_2m%i_beta_%i_eta_%.3f.txt" pdir func gamma[g] D β η
        npath2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i_eta_%.3f.txt" dir gamma[g] D β η
        mv(path2, npath2)
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