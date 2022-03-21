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
    "--width"
        arg_type = Int64
        default = 1
        help = "width of the cylinder, should be an odd number"
    "--beta"
        arg_type = Float64
        default = 1.0
        help = "β: inverse temperature"
    "--init"
        arg_type = Float64
        default = 0.0
        help = "initiate with the cMPS at β = init"            
    "--bondD"
        arg_type = Int64
        default = 8
        help = "cMPS bond dimension"
    "--power_step"
        arg_type = Int64
        default = 100
        help = "maximum power step if power method is used"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/JuliaCMPO/XXZ_2D_helical"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const Jxy = parsed_args[:Jxy]
const Jz = parsed_args[:Jz]
const width = parsed_args[:width]
const β = parsed_args[:beta]
const bondD = parsed_args[:bondD]
const power_step = parsed_args[:power_step]


#Creat ResultFolders if there is none
ResultFolder = parsed_args[:ResultFolder]

ResultFolder = "/data/sliang/JuliaCMPO/XXZ_2D_helical"
isdir(ResultFolder) || mkdir(ResultFolder)

ModelResultFolder = @sprintf "%s/Jxy_%.2f_Jz_%.2f_wid_%02i" ResultFolder Jxy Jz width
isdir(ModelResultFolder) || mkdir(ModelResultFolder)

#CMPO
model = XXZmodel_2D_helical(Jz/Jxy, width)


@timeit to "evaluate" begin
    ψ, dict = evaluate(model, bondD, β, ModelResultFolder, hermitian = false, max_pow_step = power_step)
end

for key in keys(dict)
    ObsvFolder = @sprintf "%s/Obsv_%s_bondD_%02i" ModelResultFolder key bondD
    isdir(ObsvFolder) || mkdir(ObsvFolder)
    ObsvFile = @sprintf "%s/beta_%.2f.txt" ObsvFolder β
    open(ObsvFile, "w") do file
        writedlm(file, [β dict[key]])
    end
end


const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
