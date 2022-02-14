using Plots
using DelimitedFiles, Printf

const DataFolder = "/data/sliang/CMPO/ising/spectrum_nnls/weightparam"
const FigureFolder = "/home/sliang/JuliaCode/mycMPO/figure_nnls"
const g = 1.0
const D = 8
β = 10.0

DataFile1 = @sprintf "%s/g_%.1f_D_%i_beta_%i_Akernel.txt" DataFolder g D β
DataFile2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i_Akernel.txt" DataFolder g D β
d1 = readdlm(DataFile1); λ1 = d1[:,1]; rn1 = d1[:,2]; sn1 = d1[:,3]
d2 = readdlm(DataFile2); λ2 = d2[:,1]; rn2 = d2[:,2]; sn2 = d2[:,3]

function add_sigle_point(num::Int64, data::AbstractArray)
    lambda_value = @sprintf "λ=%.5e" data[num, 1]
    scatter!([data[num, 2]], [data[num,3]]
        ,marker=(:circle, 5, :white, stroke(1.0, :red)), label=lambda_value)
end

plot(rn1, sn1, line=(2),label="D=8"
    ,xaxis=:log, xlabel="Residual(K*A-y) norm"
    ,yaxis=:log, ylabel="Solution(A) norm"
    ,framestyle=:box)
add_sigle_point(55, d1)
add_sigle_point(100, d1)
add_sigle_point(130, d1)
add_sigle_point(150, d1)
FigureName = @sprintf "ising_Lcurve_g_%.1f_D_%i_beta_%i_Akernel.pdf" g D β
savefig(FigureFolder*"/"*FigureName)


const LiFolder = "/data/sliang/CMPO/ising/spectrum_Li"
const NNLSFolder = "/data/sliang/CMPO/ising/spectrum_nnls/"
func = "Sw"
LiPath = @sprintf "%s/%s/g_%.1f_beta_%i.txt" LiFolder func g β
sli = readdlm(LiPath)

function add_line(num::Int64, data::AbstractArray)
    lambda_value = @sprintf "λ=%.5e" data[num, 1]
    path = @sprintf "%s/%s/g_%.1f_D_%i_beta_%i_lambda_%.5e_Akernel.txt" NNLSFolder func g D β data[num, 1]
    d = readdlm(path)
    plot!(d[:,1], d[:,2], label=lambda_value
            ,line=(2))
end

plot(sli[:,1], sli[:,2], label="actual"
        ,xlabel="ω", ylabel="S(ω)", xlim=(0, 5)
       # ,ylim = (0,10)
        ,line=(:black, 2))
add_line(55, d1)
add_line(100, d1)
add_line(130,d1)
add_line(150, d1)
FigureName = @sprintf "ising_Sw_g_%.1f_D_%i_beta_%i_Akernel.pdf" g D β
savefig(FigureFolder*"/"*FigureName)



plot(rn2, sn2, line=(2),label="D=8×2"
    ,xaxis=:log, xlabel="Residual(K*S-y) norm"
    ,yaxis=:log, ylabel="Solution(S) norm")