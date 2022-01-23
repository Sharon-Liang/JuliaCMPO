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
        default = "/data/sliang/CMPO/ising/spectrum_nnls/weightparam"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args)

const g = parsed_args[:g]
const β = parsed_args[:beta]
const D = parsed_args[:D]
const λ = parsed_args[:lambda]
const dω = parsed_args[:dw]
const ωmax = parsed_args[:wmax]
const DataFolder = parsed_args[:DataFolder]
const ResultFolder = parsed_args[:ResultFolder]

DataFile = @sprintf "%s/gtau/g_%.1f_D_%i_beta_%i.txt" DataFolder g D β
ResultFile = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" ResultFolder g D β

#Existence of DataFile
@assert isfile(DataFile) 

#Remove old ResultFile in the first step
if λ == 0. isfile(ResultFile) && rm(ResultFile) end

#LOAD DATA
d1 = readdlm(DataFile)
x1 = d1[:,1]; y1 = d1[:,2]

#BUILD KERNEL
ω, len = build_range(dω, ωmax)
K0 = build_kernal(0, x1, ω, β)

#NNT METHOD 
eye = Matrix(1.0I, len, len)
ȳ1 = vcat(y1, zeros(len))

K̄ = vcat(K0, λ*eye)
S = nonneg_lsq(K̄,ȳ1;alg=:nnls) 
ȳ = K0 * S

n1 = norm(ȳ)
n2 = norm(S)

#WRITE RESULT 
open(ResultFile, "a") do file
    writedlm(file, [λ n1 n2])
end

    
