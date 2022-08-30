using LinearAlgebra; BLAS.set_num_threads(Threads.nthreads())
using TimerOutputs, Dates, ArgParse
using DelimitedFiles, HDF5, Printf
using JuliaCMPO

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
    "--bi"
        arg_type = Float64
        default = 1.0
        help = "initial value of β: inverse temperature"
    "--expand"
        arg_type = Bool
        default = false
        help = "expand a cmpo or not"
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
    "--init"
        arg_type = Float64
        default = 0.
        help ="init the cmps with the one at β=init"
    "--Continue"
        arg_type = Bool
        default = false
        help ="Continue"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/JuliaCMPO/XXZ"
        help = "result folder"
    "--tag"
        arg_type = String
        default = Dates.format(now(), "yyyy-mm-dd")
        help = "date tag"
    "--processor"
        arg_type = Processor
        default = CPU
        help = "processor used"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const Jxy = parsed_args[:Jxy]
const Jz = parsed_args[:Jz]
const bi = parsed_args[:bi]
const bf = parsed_args[:bf]
const bstep = parsed_args[:bstep]
const β0 = parsed_args[:init]
const bondD = parsed_args[:bondD]
const init = parsed_args[:init]
const Continue = parsed_args[:Continue]
const tag = parsed_args[:tag]
const processor = parsed_args[:processor]

const wid = 1
const ResultFolder = parsed_args[:ResultFolder]
isdir(ResultFolder) || mkdir(ResultFolder)
ModelResultFolder = @sprintf "%s/Jz_%.2f_Jxy_%.2f_wid_%02i" ResultFolder Jz Jxy wid
isdir(ModelResultFolder) || mkdir(ModelResultFolder)


#CMPO
model = XXZmodel(Jz/Jxy)

if β0 == 0
    ψ0 = nothing
else
    init_path = @sprintf "%s/bondD_%02i_CMPS/beta_%.2f.hdf5" ModelResultFolder bondD β0
    ψ0 = readCMPS(init_path)
end

EngFile = @sprintf "%s/Obsv_FECvS_bondD_%02i.txt" ModelResultFolder bondD
if Continue == false
    open(EngFile,"w") do file  
        #write(file, "step    free_energy/site        energy/site         entropy/site    \n")
        #write(file, "----  -------------------  --------------------  -------------------\n")
        write(file, "step    free_energy/site \n")
        write(file, "----  -------------------\n")
    end
end

βlist = [i for i in range(bi, bf, step=bstep)]

for b = 1:length(βlist)
    β = βlist[b]
    @timeit to "evaluate" begin
        evaluate_options = EvaluateOptions(EvaluateOptions(),
                            init = ψ0,
                            group = group,
                            tag = tag,
                            processor = processor)
        res = JuliaCMPO.evaluate(model, bondD, β, ModelResultFolder; evaluate_options)
    end

    open(EngFile,"a") do file  
        #EngString = @sprintf "%.2f   %.16f   %.16f   %.16f \n" β res[2]["F"] res[2]["E"] res[2]["S"]
        EngString = @sprintf "%.2f   %.16f \n" β res[2]["F"]
        write(file, EngString)
    end

    global ψ0 = res[1]
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")
