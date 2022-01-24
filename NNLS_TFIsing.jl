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
        default = "/data/sliang/CMPO/ising/spectrum_nnls/weightparam"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const g = parsed_args[:g]
const β = parsed_args[:beta]
const D = parsed_args[:D]
const dω = parsed_args[:dw]
const ωmax = parsed_args[:wmax]
const DataFolder = parsed_args[:DataFolder]
const ResultFolder = parsed_args[:ResultFolder]

DataFile1 = @sprintf "%s/gtau/g_%.1f_D_%i_beta_%i.txt" DataFolder g D β
DataFile2 = @sprintf "%s/gtau/g_%.1f_D_%im2_beta_%i.txt" DataFolder g D β
ResultFile1 = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" ResultFolder g D β
ResultFile2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i.txt" ResultFolder g D β

#Existence of DataFile
@assert isfile(DataFile1) && isfile(DataFile2)

#Remove old ResultFile in the first step
isfile(ResultFile1) && rm(ResultFile1) 
isfile(ResultFile2) && rm(ResultFile2)

#LOAD DATA
d1 = readdlm(DataFile1); x1 = d1[:,1]; y1 = d1[:,2]
d2 = readdlm(DataFile2); x2 = d2[:,1]; y2 = d2[:,2]

#BUILD KERNEL
ω, len = build_range(dω, ωmax)
K1 = build_kernal(0, x1, ω, β)
K2 = build_kernal(0, x2, ω, β)

#NNT METHOD 
eye = Matrix(1.0I, len, len)
L = 300
lambda = [10^i for i in range(-5,1,length=L)]
n1, n2 = zeros(L,2), zeros(L,2)
for i=1:L
    λ = lambda[i]
    ȳ1 = vcat(y1, zeros(len)); K̄1 = vcat(K1, λ*eye)
    ȳ2 = vcat(y2, zeros(len)); K̄2 = vcat(K2, λ*eye)

    @timeit to "D=8" begin
        S1 = nonneg_lsq(K̄1,ȳ1;alg=:nnls)  
    end

    @timeit to "D=8×2" begin
        S2 = nonneg_lsq(K̄2,ȳ2;alg=:nnls) 
    end

    n1[i,1] = norm(K1 * S1 - y1); n2[i,1] = norm(S1)
    n1[i,2] = norm(K2 * S2 - y2); n2[i,2] = norm(S2)
end

#WRITE RESULT 
open(ResultFile1, "w+") do file
    for i = 1:L
        writedlm(file, [lambda[i] n1[i,1] n2[i,1]])
    end
end

open(ResultFile2, "w+") do file
    for i = 1:L
        writedlm(file, [lambda[i] n1[i,2] n2[i,2]])
    end
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
