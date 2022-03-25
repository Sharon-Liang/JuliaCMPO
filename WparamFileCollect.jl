using DelimitedFiles, Printf

const g = 1.0
const D = 8
const DataFolder = "/data/sliang/CMPO/ising/spectrum_nnls/weightparam"

lambdaRange = [10^i for i in range(-20,0,length=200)]
#betaRange = [10.0, 20.0, 30.0, 40.0]
betaRange = [10.0]

#Check the Existence of DataFiles
for b in 1:length(betaRange)
    β = betaRange[b]
    DataCollect1 = zeros(length(lambdaRange), 3)
    DataCollect2 = zeros(length(lambdaRange), 3)
    ResultFile1 = @sprintf "%s/g_%.1f_D_%i_beta_%i_dG.txt" DataFolder g D β 
    #ResultFile2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i_dG.txt" DataFolder g D β
    
    for i in 1:length(lambdaRange)
        λ = lambdaRange[i]
        DataFile1 = @sprintf "%s/g_%.1f_D_%i_beta_%i_lambda_%.5e_dG.txt" DataFolder g D β λ
        DataFile2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i_lambda_%.5e_dG.txt" DataFolder g D β λ
        DataCollect1[i, :] = readdlm(DataFile1)
        #DataCollect2[i, :] = readdlm(DataFile2)
    end

    open(ResultFile1, "w+") do file
        writedlm(file, DataCollect1)
    end
    #open(ResultFile2, "w+") do file
    #    writedlm(file, DataCollect2)
    #end
end



#CLEAR LOG FOLDER
#if length(readdir(logdir))!=0
#    for file in readdir(logdir)
#    run(```rm $(logdir)/$(file)```) end
#end

using Printf
ResultFolder = "/data/sliang/CMPO/TFIsing/Power_J_1.00_G_1.00/"
χ = 8
for β in range(1.0, 16.7, step = 0.1)
    #β = 16.9
    OldFile = @sprintf "%s/bondD_%2i_beta_%.1f.hdf5" ResultFolder χ β
    NewFile = @sprintf "%s/bondD_%02i_beta_%.2f.hdf5" ResultFolder χ β
    #run(```cp $(OldFile) $(NewFile)```)
    run(```rm $(OldFile)```)
end


using Printf

for wid in [1,3,5]
    for bondD in [8, 10, 12]
        for β in [1, 2, 4, 8, 16]
            ChkpFolder = @sprintf "/data/sliang/JuliaCMPO/XXZ_2D_helical/Jz_1.00_Jxy_1.00_wid_%02i/bondD_%02i_CMPS/CheckPoint_beta_%.2f" wid bondD β

            if isdir(ChkpFolder)
                println(ChkpFolder)
                FFile = @sprintf "%s/Obsv_F_fidelity.txt" ChkpFolder
                if isfile(FFile)
                    rm(FFile)
                end
            end
        end
    end
end