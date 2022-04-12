using LinearAlgebra; BLAS.set_num_threads(Threads.nthreads())
using TimerOutputs, Dates, ArgParse
using DelimitedFiles, HDF5, Printf
using cMPO

println("JULIA_NUM_THREADS = ", Threads.nthreads())
println("OPENBLAS_NUM_THREADS = ", BLAS.get_num_threads())

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
    "--expand"
        arg_type = Bool
        default = false
        help = "expand a cmpo or not"
    "--group"
        arg_type = Int64
        default = 2
        help = "group :group sites as a unit cell"
    "--beta"
        arg_type = Float64
        default = 1.0
        help = "β: inverse temperature"          
    "--bondD"
        arg_type = Int64
        default = 8
        help = "cMPS bond dimension"
    "--max_pow_step"
        arg_type = Int64
        default = 100
        help = "maximum power step if power method is used"
    "--Continue"
        arg_type = Int64
        default = 0
        help ="Continue = 0 :false, Continue > max_pow_step : true, else Continue = Continue"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/JuliaCMPO/XXZ_2D_helical"
        help = "result folder"
    "--tag"
        arg_type = String
        default = Dates.format(now(), "yyyy-mm-dd")
        help = "date tag"
    "--device"
        arg_type = Symbol
        default = :cpu
        help = "device used"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const Jxy = parsed_args[:Jxy]
const Jz = parsed_args[:Jz]
const width = parsed_args[:width]
const expand = parsed_args[:expand]
const group = parsed_args[:group]
const β = parsed_args[:beta]
const bondD = parsed_args[:bondD]
const max_pow_step = parsed_args[:max_pow_step]
const tag = parsed_args[:tag]
const device = parsed_args[:device]

Continue = parsed_args[:Continue]
#Creat ResultFolders if there is none
const ResultFolder = parsed_args[:ResultFolder]
isdir(ResultFolder) || mkdir(ResultFolder)

ModelResultFolder = @sprintf "%s/Jz_%.2f_Jxy_%.2f_wid_%02i" ResultFolder Jz Jxy width
isdir(ModelResultFolder) || mkdir(ModelResultFolder)

#CMPO
model = XXZmodel_2D_helical(Jz/Jxy, width, expand = expand)

if Continue ≥ max_pow_step Continue = true end
@timeit to "evaluate" begin
    res = evaluate(model, bondD, β, ModelResultFolder, 
                        hermitian = false, 
                        max_pow_step = max_pow_step, 
                        group = group,
                        Continue = Continue,
                        tag = tag,
                        device=device)
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
