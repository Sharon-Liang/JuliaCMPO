using TimerOutputs, Dates, ArgParse
using DelimitedFiles, Printf
using NonNegLeastSquares, LinearAlgebra
include("/home/sliang/JuliaCode/mycMPO/src/KernelFunction.jl")

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="NNLS code for XXZ model"
)
@add_arg_table! settings begin
    "--Jxy"
        arg_type = Float64
        default = 1.0
        help = "Jxy: coupling constant for Sx and Sy terms"
    "--Jz"
        arg_type = Float64
        default = 1.0
        help = "Jz: coupling constant for Sz term"
    "--beta"
        arg_type = Float64
        default = 10.0
        help = "β: inverse temperature"
    "--D"
        arg_type = Int64
        default = 8
        help = "cMPO bond dimension"
    "--operator"
        arg_type = Symbol
        default = :mp
        help = "spin operators in the correlation function: :mp, :pm, :px, :py, :pz"
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
        default = "/data/sliang/CMPO/xxz/imagtime"
        help = "data folder"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/CMPO/xxz/spectrum_nnls"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const Jxy = parsed_args[:Jxy]
const Jz = parsed_args[:Jz]
const β = parsed_args[:beta]
const D = parsed_args[:D]
const op = parsed_args[:operator]
const λ = parsed_args[:lambda]
const dω = parsed_args[:dw]
const ωmax = parsed_args[:wmax]
const DataFolder = parsed_args[:DataFolder]
const ResultFolder = parsed_args[:ResultFolder]

#DataFile Path 
DataFile1 = @sprintf "%s/gtau/Jz_%.1f_%s_D_%i_beta_%i.txt" DataFolder Jz op D β
DataFile2 = @sprintf "%s/gtau/Jz_%.1f_%s_D_%im2_beta_%i.txt" DataFolder Jz op D β
#Existence of DataFile
println(DataFile1)
@assert isfile(DataFile1) && isfile(DataFile2)

#Creat ResultFolders if there is none
isdir(ResultFolder) || mkdir(ResultFolder)
isdir(ResultFolder*"/weightparam") || mkdir(ResultFolder*"/weightparam")
isdir(ResultFolder*"/Sw") || mkdir(ResultFolder*"/Sw")
isdir(ResultFolder*"/Aw") || mkdir(ResultFolder*"/Aw")

#ResultFile Path
WP_ResultFile1 = @sprintf "%s/weightparam/Jz_%.1f_%s_D_%i_beta_%i_lambda_%e.txt" ResultFolder Jz op D β λ
WP_ResultFile2 = @sprintf "%s/weightparam/Jz_%.1f_%s_D_%im2_beta_%i_lambda_%e.txt" ResultFolder Jz op D β λ

SW_ResultFile1 = @sprintf "%s/Sw/Jz_%.1f_%s_D_%i_beta_%i_lambda_%e.txt" ResultFolder Jz op D β λ
SW_ResultFile2 = @sprintf "%s/Sw/Jz_%.1f_%s_D_%im2_beta_%i_lambda_%e.txt" ResultFolder Jz op D β λ

AW_ResultFile1 = @sprintf "%s/Aw/Jz_%.1f_%s_D_%i_beta_%i_lambda_%e.txt" ResultFolder Jz op D β λ
AW_ResultFile2 = @sprintf "%s/Aw/Jz_%.1f_%s_D_%im2_beta_%i_lambda_%e.txt" ResultFolder Jz op D β λ

#Remove old ResultFile in the first step
isfile(WP_ResultFile1) && rm(WP_ResultFile1) 
isfile(WP_ResultFile2) && rm(WP_ResultFile2)

isfile(SW_ResultFile1) && rm(SW_ResultFile1) 
isfile(SW_ResultFile2) && rm(SW_ResultFile2)

isfile(AW_ResultFile1) && rm(AW_ResultFile1) 
isfile(AW_ResultFile2) && rm(AW_ResultFile2)


#LOAD DATA
d1 = readdlm(DataFile1); x1 = d1[:,1]; y1 = d1[:,2]
d2 = readdlm(DataFile2); x2 = d2[:,1]; y2 = d2[:,2]

#BUILD KERNEL
ω, len = build_range(dω, ωmax)
K1 = build_kernal(0, x1, ω, β)
K2 = build_kernal(0, x2, ω, β)

#NNT METHOD 
eye = Matrix(1.0I, len, len)

ȳ1 = vcat(y1, zeros(len)); K̄1 = vcat(K1, λ*eye)
ȳ2 = vcat(y2, zeros(len)); K̄2 = vcat(K2, λ*eye)

@timeit to "D=8" begin
    S1 = nonneg_lsq(K̄1,ȳ1;alg=:nnls)
    A1 = map((x,y) -> ASrelation(x,y,β), ω, S1)
end

@timeit to "D=8×2" begin
    S2 = nonneg_lsq(K̄2,ȳ2;alg=:nnls) 
    A2 = map((x,y) -> ASrelation(x,y,β), ω, S2)
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

open(SW_ResultFile1, "w+") do file
    for i = 1:len
        writedlm(file, [ω[i] S1[i]])
    end
end

open(SW_ResultFile2, "w+") do file
    for i = 1:len
        writedlm(file, [ω[i] S2[i]])
    end
end

open(AW_ResultFile1, "w+") do file
    for i = 1:len
        writedlm(file, [ω[i] A1[i]])
    end
end

open(AW_ResultFile2, "w+") do file
    for i = 1:len
        writedlm(file, [ω[i] A2[i]])
    end
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
