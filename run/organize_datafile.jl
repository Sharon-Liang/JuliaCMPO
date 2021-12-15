using DelimitedFiles
using JLD, HDF5
using Printf


# rename file
model = "ising"
beta =[1.0, 2.0, 4.0, 6.0, 8.0, 10., 20., 30., 40.]
gamma=[1.0]


folder = "spectrum-Li"
func = "S0iwn"
pdir = @sprintf "./data/%s/%s/data/change_eta" model folder
dir = @sprintf "./data/%s/%s/%s" model folder func
isdir(dir) || mkdir(dir)

for g = 1:length(gamma)
    (d0, r0) = divrem(g, 1)
    for b = 1:length(beta)
        β = beta[b]
        path1 = @sprintf "%s/g_%ip%i_beta_%i.txt" pdir d0 r0 β
        npath1 = @sprintf "%s/g_%.1f_beta_%i.txt" dir gamma[g] β 
        mv(path1, npath1)
        #path2 = @sprintf "%s/g_%ip%i_beta_%i_%s_eta_2pidbeta.txt" pdir d0 r0 β func
        #npath2 = @sprintf "%s/g_%.1f_beta_%i_eta_2pidbeta.txt" dir gamma[g] β 
        #mv(path2, npath2)
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