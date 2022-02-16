using TimerOutputs, Dates, ArgParse
using DelimitedFiles, JLD, HDF5, Printf
using cMPO
using Optim

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="NNLS code for TFIsing model"
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
    "--beta"
        arg_type = Float64
        default = 10.0
        help = "β: inverse temperature"
    "--prevbeta"
        arg_type = Float64
        default = 0.0
        help = "using the result wave fuction of β = prevbeta as initial wave function"
    "--D"
        arg_type = Int64
        default = 8
        help = "cMPO bond dimension"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/CMPO/ising"
        help = "result folder"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const J = parsed_args[:J]
const Γ = parsed_args[:gamma]
const β = parsed_args[:beta]
const prevbeta = parsed_args[:prevbeta]
const D = parsed_args[:D]
const ResultFolder = parsed_args[:ResultFolder]

#Creat ResultFolders if there is none
isdir(ResultFolder) || mkdir(ResultFolder)

#ResultFile Path
if J == 1.0
    ResultFile = @sprintf "%s/g_%.1f_D_%i_beta_%i.jld" ResultFolder Γ D β
    DataFile = @sprintf "%s/g_%.1f_D_%i_beta_%i.jld" ResultFolder Γ D prevbeta
else
    ResultFile = @sprintf "%s/j_%.1f_D_%i_beta_%i.jld" ResultFolder J D β
    DataFile = @sprintf "%s/j_%.1f_D_%i_beta_%i.jld" ResultFolder Γ D prevbeta
end

#Remove old ResultFile in the first step
isfile(ResultFile) && rm(ResultFile) 

#CMPO
w = TFIsing(J, Γ)

if prevbeta == 0.0
    arr = init_cmps(D, w) |> toarray 
else
    arr = load(DataFile)[string(prevbeta)][2]
end
f = arr -> free_energy(arr, w, β)
gf! = gradient_function(f)
opt = optimize(f,gf!, arr, LBFGS(),
        Optim.Options(iterations = 10000, store_trace=true, time_limit = 7200))

      
#WRITE RESULT
arr = opt.minimizer
key = string(β)
if Optim.converged(opt) == true
    res = (minimum(opt), arr, Optim.converged(opt))
else
    res = (minimum(opt), arr, Optim.converged(opt), Optim.f_trace(opt))
end

jldopen(Res, "w+") do file
    write(ResultFile, key, res)
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")

