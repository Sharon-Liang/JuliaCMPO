using LinearAlgebra; BLAS.set_num_threads(Threads.nthreads())
using TimerOutputs, Dates, ArgParse
using DelimitedFiles, Printf
using JuliaCMPO

println("JULIA_NUM_THREADS = ", Threads.nthreads())
println("OPENBLAS_NUM_THREADS = ", BLAS.get_num_threads())

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="CMPO: 1D TFIsing model"
)
@add_arg_table! settings begin
    "--J1"
        arg_type = Float64
        default = 1.0
        help = "Transverse Field Ising σz coupling constant"
    "--J2"
        arg_type = Float64
        default = 1.0
        help = "Transverse Field strength: Γ"
    "--beta"
        arg_type = Float64
        default = 1.0
        help = "β: inverse temperature"
    "--bondD"
        arg_type = Int64
        default = 8
        help = "cMPS bond dimension"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/JuliaCMPO/J1J2Chain"
        help = "result folder"
    "--processor"
        arg_type = Int64
        default = 0
        help = "processor used: 0 for CPU and 1 for GPU"
    "--tag"
        arg_type = String
        default = Dates.format(now(), "yyyy-mm-dd")
        help = "date tag"
    "--group"
        arg_type = Int
        default = 0
        help = "size of unit cell"
    "--shift"
        arg_type = Float64
        default = 0.
        help = "shift spectrum"
end

parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const J1 = parsed_args[:J1]
const J2 = parsed_args[:J2]
const β = parsed_args[:beta]
const bondD = parsed_args[:bondD]
const ResultFolder = parsed_args[:ResultFolder]
const tag = parsed_args[:tag]
const processor = Processor(parsed_args[:processor])

const to_group = parsed_args[:group]
const to_shift = parsed_args[:shift]


const wid = 1

isdir(ResultFolder) || mkdir(ResultFolder)

ModelResultFolder = @sprintf "%s/J1_%f_J2_%f_wid_%02i_bondD_%02i_group_%i_shift_%e_%s" ResultFolder J1 J2 wid bondD to_group to_shift tag 
isdir(ModelResultFolder) || mkdir(ModelResultFolder)

ModelResultFolder = @sprintf "%s/beta_%f" ModelResultFolder β
isdir(ModelResultFolder) || mkdir(ModelResultFolder)
   
#CMPO
Tₘ = model(J1J2Chain(), J1, J2)

@timeit to "evaluate" begin
    power_evaluate(Tₘ, bondD, [β]; processor, result_folder = ModelResultFolder, to_shift, to_group)
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)

@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")

