using LinearAlgebra; BLAS.set_num_threads(Threads.nthreads())
using TimerOutputs, Dates, ArgParse
using DelimitedFiles, HDF5, Printf
using JuliaCMPO

println("JULIA_NUM_THREADS = ", Threads.nthreads())
println("OPENBLAS_NUM_THREADS = ", BLAS.get_num_threads())

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="CMPO: TFIsing model"
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
        default = "/data/sliang/JuliaCMPO/TFIsing"
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

const J = parsed_args[:J]
const Γ = parsed_args[:G]
const bi = parsed_args[:bi]
const bf = parsed_args[:bf]
const bstep = parsed_args[:bstep]
const β0 = parsed_args[:init]
const bondD = parsed_args[:bondD]
const Continue = parsed_args[:Continue]
const ResultFolder = parsed_args[:ResultFolder]
const tag = parsed_args[:tag]
const processor = parsed_args[:processor]

const wid = 1

isdir(ResultFolder) || mkdir(ResultFolder)
ModelResultFolder = @sprintf "%s/J_%.2f_G_%.2f_wid_%02i" ResultFolder J Γ wid
isdir(ModelResultFolder) || mkdir(ModelResultFolder)
   
#CMPO
model = TFIsing(J, Γ)

if β0 == 0
    ψ0 = nothing
else
    init_path = @sprintf "%s/bondD_%02i_CMPS/beta_%.2f.hdf5" ModelResultFolder bondD β0
    ψ0 = readCMPS(init_path)
end

EngFile = @sprintf "%s/Obsv_FECvS_bondD_%02i.txt" ModelResultFolder bondD
if Continue == false
    open(EngFile,"w") do file  
        write(file, "  β          free_energy           energy              specific_heat            entropy      \n")
        write(file, "------  -------------------  --------------------   -------------------  -------------------\n")
    end
end

βlist = [i for i in range(bi, bf, step=bstep)]

for b = 1:length(βlist)
    β = βlist[b]
    @timeit to "evaluate" begin
        evaluate_options = EvaluateOptions(EvaluateOptions(),
                            init = ψ0,
                            tag = tag,
                            processor = processor)
        res = JuliaCMPO.evaluate(model, bondD, β, ModelResultFolder, 
                       options = evaluate_options)
    end

    open(EngFile,"a") do file  
        EngString = @sprintf "%.2f   %.16f   %.16f   %.16f   %.16f \n" β res[2]["F"] res[2]["E"] res[2]["Cv"] res[2]["S"]
        write(file, EngString)
    end

    global ψ0 = res[1]
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")

