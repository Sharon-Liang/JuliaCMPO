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
        arg_type = Int64
        default = 8
        help = "cMPO bond dimension"
    "--lambda"
        arg_type = Float64
        default = 0.0
        help = "weight parameter"
    "--spectrum"
        arg_type = Symbol
        default = :S
        help = "spectrum to be calculated"
    "--dw"
        arg_type = Float64
        default = 0.01
        help = "ω step size"
    "--wmin"
        arg_type = Float64
        default = 0.0
        help = "minimum ω value"
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
const spectrum = parsed_args[:spectrum]
const dω = parsed_args[:dw]
const ωmin = parsed_args[:wmin]
const ωmax = parsed_args[:wmax]
const DataFolder = parsed_args[:DataFolder]
const ResultFolder = parsed_args[:ResultFolder]

#DataFile Path 
DataFile1 = @sprintf "%s/gtau/g_%.1f_D_%i_beta_%i.txt" DataFolder g D β
DataFile2 = @sprintf "%s/gtau/g_%.1f_D_%im2_beta_%i.txt" DataFolder g D β
#Existence of DataFile
@assert isfile(DataFile1) && isfile(DataFile2)

#Creat ResultFolders if there is none
isdir(ResultFolder) || mkdir(ResultFolder)
isdir(ResultFolder*"/weightparam") || mkdir(ResultFolder*"/weightparam")
isdir(ResultFolder*"/weightparam/$(spectrum)w") || mkdir(ResultFolder*"/weightparam/$(spectrum)w")
isdir(ResultFolder*"/$(spectrum)w") || mkdir(ResultFolder*"/$(spectrum)w")

#ResultFile Path
WP_ResultFile1 = @sprintf "%s/weightparam/%sw/g_%.1f_D_%i_beta_%i_lambda_%e.txt" ResultFolder spectrum g D β λ
WP_ResultFile2 = @sprintf "%s/weightparam/%sw/g_%.1f_D_%im2_beta_%i_lambda_%e.txt" ResultFolder spectrum g D β λ

ResultFile1 = @sprintf "%s/%sw/g_%.1f_D_%i_beta_%i_lambda_%e.txt" ResultFolder spectrum g D β λ
ResultFile2 = @sprintf "%s/%sw/g_%.1f_D_%im2_beta_%i_lambda_%e.txt" ResultFolder spectrum g D β λ

#Remove old ResultFile in the first step
isfile(WP_ResultFile1) && rm(WP_ResultFile1) 
isfile(WP_ResultFile2) && rm(WP_ResultFile2)

isfile(ResultFile1) && rm(ResultFile1) 
isfile(ResultFile2) && rm(ResultFile2)


#LOAD DATA
d1 = readdlm(DataFile1); x1 = d1[:,1]; y1 = d1[:,2]
d2 = readdlm(DataFile2); x2 = d2[:,1]; y2 = d2[:,2]

#BUILD KERNEL
ω, len = build_range(dω, ωmax, ωmin)
K1 = build_kernal(0, x1, ω, β, spectrum)
K2 = build_kernal(0, x2, ω, β, spectrum)

#NNT METHOD 
eye = Matrix(1.0I, len, len)

ȳ1 = vcat(y1, zeros(len)); K̄1 = vcat(K1, λ*eye)
ȳ2 = vcat(y2, zeros(len)); K̄2 = vcat(K2, λ*eye)

@timeit to "D=8" begin
    S1 = nonneg_lsq(K̄1,ȳ1;alg=:nnls)
end

@timeit to "D=8×2" begin
    S2 = nonneg_lsq(K̄2,ȳ2;alg=:nnls) 
end

Rnorm1 = norm(K1 * S1 - y1); Snorm1 = norm(S1)
Rnorm2 = norm(K2 * S2 - y2); Snorm2 = norm(S2)


#WRITE RESULT 
open(WP_ResultFile1, "w+") do file
    writedlm(file, [λ Rnorm1 Snorm1])
end

open(WP_ResultFile2, "w+") do file
    writedlm(file, [λ Rnorm2 Snorm2])
end

open(ResultFile1, "w+") do file
    for i = 1:len
        writedlm(file, [ω[i] S1[i]])
    end
end

open(ResultFile2, "w+") do file
    for i = 1:len
        writedlm(file, [ω[i] S2[i]])
    end
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
