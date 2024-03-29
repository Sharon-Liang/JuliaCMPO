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
    "--J"
        arg_type = Float64
        default = 1.0
        help = "Transverse Field Ising σz coupling constant"
    "--G"
        arg_type = Float64
        default = 1.0
        help = "Transverse Field strength: Γ"
    "--bi"
        arg_type = Float64
        default = 1.0
        help = "initial value of β: inverse temperature"
    "--bf"
        arg_type = Float64
        default = 30.0
        help = "final value of β: inverse temperature" 
    "--bstep"
        arg_type = Float64
        default = 1.0
        help = "step size of β: inverse temperature" 
    "--bondD"
        arg_type = Int64
        default = 8
        help = "cMPS bond dimension"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/JuliaCMPO/TFIsing"
        help = "result folder"
    "--processor"
        arg_type = Int64
        default = 0
        help = "processor used: 0 for CPU and 1 for GPU"
    "--tag"
        arg_type = String
        default = Dates.format(now(), "yyyy-mm-dd")
        help = "date tag"
end

parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const J = parsed_args[:J]
const Γ = parsed_args[:G]
const bi = parsed_args[:bi]
const bf = parsed_args[:bf]
const bstep = parsed_args[:bstep]
const bondD = parsed_args[:bondD]
const ResultFolder = parsed_args[:ResultFolder]
const tag = parsed_args[:tag]
const processor = Processor(parsed_args[:processor])


const wid = 1

isdir(ResultFolder) || mkdir(ResultFolder)
ModelResultFolder = @sprintf "%s/J_%.2f_G_%.2f_wid_%02i_bondD_%02i_%s" ResultFolder J Γ wid bondD tag
isdir(ModelResultFolder) || mkdir(ModelResultFolder)
   
#CMPO
Tₘ = model(TFIsingChain(), Γ, J)
βlist = [i for i in range(bi, bf, step=bstep)]

@timeit to "evaluate" begin
    variation_evaluate(Tₘ, bondD, βlist; processor, result_folder = ModelResultFolder)
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)

@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")

