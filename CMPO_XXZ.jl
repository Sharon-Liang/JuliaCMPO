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
    "--D"
        arg_type = Int64
        default = 8
        help = "cMPO bond dimension"
    "--width"
        arg_type = Int64
        default = 1
        help = "width of the cylinder, width = 1 for 1D chain"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/CMPO/xxz"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const Jxy = parsed_args[:Jxy]
const Jz = parsed_args[:Jz]
const beta_max = parsed_args[:beta_max]
const beta_min = parsed_args[:beta_min]
const beta_step = parsed_args[:beta_step]
const D = parsed_args[:D]
const width = parsed_args[:width]
const ResultFolder = parsed_args[:ResultFolder]

#Creat ResultFolders if there is none
isdir(ResultFolder) || mkdir(ResultFolder)



#ResultFile Path
isdir(ResultFolder*"/FreeEnergy") || mkdir(ResultFolder*"/FreeEnergy")
Eng_ResultFile = @sprintf "%s/FreeEnergy/Jz_%.1f_W_%i_D_%i.txt" ResultFolder Jz width D
#Remove old ResultFile in the first step
isfile(Eng_ResultFile) && rm(Eng_ResultFile) 

#CMPO
Δ = Jz/Jxy
w = XXZmodel(Δ, width)
beta_range = [i for i in range(beta_min, beta_max, step = beta_step)]
Fenergy = similar(beta_range)

arr = init_cmps(D, w) |> toarray
@show free_energy(arr, w, beta_min)
for b = 1:length(beta_range)
    β = beta_range[b]; key = string(β)  
    ResultFile = @sprintf "%s/Jz_%.1f_W_%i_D_%i_beta_%.2f.jld" ResultFolder Jz width D β
    isfile(ResultFile) && rm(ResultFile) 

    f = x -> free_energy(x, w, β)
    gf! = gradient_function(f)
    @timeit to "β = $(key)" begin
        opt = optimize(f,gf!, arr, LBFGS(),
                Optim.Options(iterations = 10000, time_limit = 7200))
    end
    #WRITE RESULT
    global arr = opt.minimizer
    Fenergy[b] = minimum(opt)
    res = (minimum(opt), arr, Optim.converged(opt))

    jldopen(ResultFile, "w") do file
        write(file, key, res)
    end

    if Optim.converged(opt) == false
        println("β = $(key)")
        @show opt
    end
end

open(Eng_ResultFile,"w+") do file
    for i=1: length(beta_range)
        writedlm(file,[beta_range[i] Fenergy[i]])
    end
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
