using TimerOutputs, Dates, ArgParse
using DelimitedFiles, JLD, HDF5, Printf
using cMPO
using Optim

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="CMPO code for XXZ model"
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
    "--beta_min"
        arg_type = Float64
        default = 1.0
        help = "minimum β(inverse temperature) value"
    "--beta_max"
        arg_type = Float64
        default = 20.0
        help = "maximum β(inverse temperature) value"
    "--beta_step"
        arg_type = Float64
        default = 0.1
        help = "β(inverse temperature) step size"
    "--init"
        arg_type = Float64
        default = 0.0
        help = "initiate with the cMPS at β = init"            
    "--bondD"
        arg_type = Int64
        default = 8
        help = "cMPS bond dimension"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/CMPO/XXZ"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const Jxy = parsed_args[:Jxy]
const Jz = parsed_args[:Jz]
const beta_max = parsed_args[:beta_max]
const beta_min = parsed_args[:beta_min]
const beta_step = parsed_args[:beta_step]
const bondD = parsed_args[:bondD]
const init = parsed_args[:init]


#Creat ResultFolders if there is none
ResultFolder = parsed_args[:ResultFolder]
isdir(ResultFolder) || mkdir(ResultFolder)

ResultFolder = @sprintf "%s/Jxy_%.2f_Jz_%.2f" ResultFolder Jxy Jz
isdir(ResultFolder) || mkdir(ResultFolder)


#CMPO
model = XXZmodel(Jz/Jxy)

if init == 0.0
    ψ0 = init_cmps(bondD, model, beta_min)
else
    initFilePath = @sprintf "%s/CMPS/bondD_%02i_beta_%.2f.hdf5" ResultFolder bondD init
    ψ0 = readCMPS(initFilePath) 
end

for β in range(beta_min, beta_max, step = beta_step)
    @timeit to "evaluate" begin
        ψ, dict = evaluate(model, bondD, β, ResultFolder, init = ψ0)
    end
    for key in keys(dict)
        filename = @sprintf "%s/Obsv_%s_bondD_%02i.txt" ResultFolder key bondD
        β == beta_min ? mode = "w" : mode = "a"
        open(filename, mode) do file
            writedlm(file, [β dict[key]])
        end
    end
    global ψ0 = ψ
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
