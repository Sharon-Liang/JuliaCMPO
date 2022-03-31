using LinearAlgebra; BLAS.set_num_threads(Threads.nthreads())
using TimerOutputs, Dates, ArgParse
using DelimitedFiles, HDF5, Printf
using cMPO

println("JULIA_NUM_THREADS = ", Threads.nthreads())
println("OPENBLAS_NUM_THREADS = ", BLAS.get_num_threads())

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="CMPO code for TFIsing_2D_helical"
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
    "--width"
        arg_type = Int64
        default = 1
        help = "width of the cylinder, should be an odd number"
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
        default = "/data/sliang/JuliaCMPO/TFIsing_2D_helical_unexpanded"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const J = parsed_args[:J]
const Γ = parsed_args[:G]
const width = parsed_args[:width]
const β = parsed_args[:beta]
const bondD = parsed_args[:bondD]
const max_pow_step = parsed_args[:max_pow_step]
Continue = parsed_args[:Continue]
#Creat ResultFolders if there is none
const ResultFolder = parsed_args[:ResultFolder]
isdir(ResultFolder) || mkdir(ResultFolder)

ModelResultFolder = @sprintf "%s/J_%.2f_G_%.2f_wid_%02i" ResultFolder J Γ width
isdir(ModelResultFolder) || mkdir(ModelResultFolder)

#CMPO
model = TFIsing_2D_helical(J, Γ, width, expand=false)
if Continue ≥ max_pow_step Continue = true end
@timeit to "evaluate" begin
    res = evaluate(model, bondD, β, ModelResultFolder, 
                        hermitian = false, 
                        max_pow_step = max_pow_step, Continue = Continue)
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
