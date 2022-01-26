using DelimitedFiles, Printf

const g = 1.0
const D = 8
const DataFolder = "/data/sliang/CMPO/ising/spectrum_nnls/weightparam"

lambdaRange = [10^i for i in range(-20,0,length=200)]
betaRange = [10.0, 20.0, 30.0, 40.0]


#Check the Existence of DataFiles
for b in 1:length(betaRange)
    β = betaRange[b]
    DataCollect1 = zeros(length(lambdaRange), 3)
    DataCollect2 = zeros(length(lambdaRange), 3)
    ResultFile1 = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" DataFolder g D β 
    ResultFile2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i.txt" DataFolder g D β
    
    for i in 1:length(lambdaRange)
        λ = lambdaRange[i]
        DataFile1 = @sprintf "%s/g_%.1f_D_%i_beta_%i_lambda_%.5e.txt" DataFolder g D β λ
        DataFile2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i_lambda_%.5e.txt" DataFolder g D β λ
        DataCollect1[i, :] = readdlm(DataFile1)
        DataCollect2[i, :] = readdlm(DataFile2)
    end

    open(ResultFile1, "w+") do file
        writedlm(file, DataCollect1)
    end
    open(ResultFile2, "w+") do file
        writedlm(file, DataCollect2)
    end
end

