using TimerOutputs, Dates, ArgParse
using DelimitedFiles, JLD, HDF5, Printf
using cMPO
using Optim

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="CMPO: TFIsing model"
)
@add_arg_table! settings begin
    "--J"
        arg_type = Float64
        default = 1.0
        help = "Transverse Field Ising σz coupling constant"
    "--gamma"
        arg_type = Float64
        default = 1.0
        help = "Transverse Field strength: Γ"
    "--beta_min"
        arg_type = Float64
        default = 1.0
        help = "minimum value of β(inverse temperature)"
    "--beta_max"
        arg_type = Float64
        default = 1.0
        help = "maximum value of β(inverse temperature)"
    "--beta_step"
        arg_type = Float64
        default = 1.0
        help = "step size of β(inverse temperature)"         
    "--bondD"
        arg_type = Int64
        default = 8
        help = "cMPS bond dimension"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/JuliaCMPO/TFIsing"
        help = "result folder"
end

parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const J = parsed_args[:J]
const Γ = parsed_args[:gamma]
const beta_min = parsed_args[:beta_min]
const beta_max = parsed_args[:beta_max]
const beta_step = parsed_args[:beta_step]
const bondD = parsed_args[:bondD]

#Creat ResultFolders if there is none
ResultFolder = parsed_args[:ResultFolder]
isdir(ResultFolder) || mkdir(ResultFolder)

ResultFolder = @sprintf "%s/Power_J_%.2f_G_%.2f_random_init" ResultFolder J Γ

   
#CMPO
model = TFIsing(J, Γ)

for β in range(beta_min, beta_max, step = beta_step)
    @timeit to "evaluate" begin
        ψ, dict = evaluate(model, bondD, β, ResultFolder, hermitian = false, init = init_cmps(bondD), max_pow_step = 50)
    end
    for key in keys(dict)
        filename = @sprintf "%s/Obsv_%s_bondD_%02i.txt" ResultFolder key bondD
        open(filename, "a") do file
            writedlm(file, [β dict[key]])
        end
    end
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")

