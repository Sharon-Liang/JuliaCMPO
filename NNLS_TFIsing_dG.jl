using TimerOutputs, Dates, ArgParse
using DelimitedFiles, Printf
using NonNegLeastSquares, LinearAlgebra
include("/home/sliang/JuliaCode/mycMPO/src/KernelFunction.jl")

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="NNLS code for TFIsing model"
)
@add_arg_table! settings begin
    "--g"
        arg_type = Float64
        default = 1.0
        help = "Transverse Field Ising Γ/J"
    "--beta"
        arg_type = Float64
        default = 10.0
        help = "β: inverse temperature"
    "--D"
        arg_type = Integer
        default = 8
        help = "cMPO bond dimension"
    "--lambda"
        arg_type = Float64
        default = 0.0
        help = "weight parameter"
    "--dw"
        arg_type = Float64
        default = 0.01
         help = "ω step size"
    "--wmax"
        arg_type = Float64
        default = 10.0
        help = "maximum ω value"
    "--DataFolder"
        arg_type = String
        default = "/data/sliang/CMPO/ising/imagtime"
        help = "data folder"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/CMPO/ising/spectrum_nnls"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const g = parsed_args[:g]
const β = parsed_args[:beta]
const D = parsed_args[:D]
const λ = parsed_args[:lambda]
const dω = parsed_args[:dw]
const ωmax = parsed_args[:wmax]
const DataFolder = parsed_args[:DataFolder]
const ResultFolder = parsed_args[:ResultFolder]

#DataFile Path 
DataFile1 = @sprintf "%s/gtau/g_%.1f_D_%i_beta_%i.txt" DataFolder g D β
DataFile2 = @sprintf "%s/dgtau/g_%.1f_D_%i_beta_%i.txt" DataFolder g D β
#DataFile4 = @sprintf "%s/dgtau/g_%.1f_D_%im2_beta_%i.txt" DataFolder g D β

#Existence of DataFile
@assert isfile(DataFile1) && isfile(DataFile2)

#ResultFile Path
WP_ResultFile1 = @sprintf "%s/weightparam/g_%.1f_D_%i_beta_%i_lambda_%e_dG.txt" ResultFolder g D β λ

SW_ResultFile1 = @sprintf "%s/Sw/g_%.1f_D_%i_beta_%i_lambda_%e_dG.txt" ResultFolder g D β λ

AW_ResultFile1 = @sprintf "%s/Aw/g_%.1f_D_%i_beta_%i_lambda_%e_dG.txt" ResultFolder g D β λ

#Remove old ResultFile in the first step
isfile(WP_ResultFile1) && rm(WP_ResultFile1) 

isfile(SW_ResultFile1) && rm(SW_ResultFile1) 

isfile(AW_ResultFile1) && rm(AW_ResultFile1) 


#LOAD DATA
d1 = readdlm(DataFile1); x1 = d1[:,1]; y1 = d1[:,2]
d2 = readdlm(DataFile2); x2 = d2[:,1]; y2 = d2[:,2]

#BUILD KERNEL
ω, len = build_range(dω, ωmax)
K1 = build_kernal(0, x1, ω, β)
K2 = build_kernal(1, x2, ω, β)


#NNT METHOD 
eye = Matrix(1.0I, len, len)

ȳ1 = vcat(y1, y2, zeros(len)); K̄1 = vcat(K1, K2, λ*eye)

@timeit to "D=8" begin
    S1 = nonneg_lsq(K̄1,ȳ1;alg=:nnls)
    A1 = map((x,y) -> ASrelation(x,y,β), ω, S1)
end

Rnorm1 = norm(vcat(K1, K2) * S1 - vcat(y1,y2)); Snorm1 = norm(S1)

#WRITE RESULT 
open(WP_ResultFile1, "w+") do file
    writedlm(file, [λ Rnorm1 Snorm1])
end

open(SW_ResultFile1, "w+") do file
    for i = 1:len
        writedlm(file, [ω[i] S1[i]])
    end
end

open(AW_ResultFile1, "w+") do file
    for i = 1:len
        writedlm(file, [ω[i] A1[i]])
    end
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
