using TimerOutputs, Dates, ArgParse
using DelimitedFiles, JLD, HDF5, Printf
using cMPO
using Optim

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="Correlations: TFIsing model"
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
        default = 1.0
        help = "β: inverse temperature"
    "--ol"
        arg_type = Symbol
        default = :z
        help = "spin operators in the correlation function <ol or>"
    "--or"
        arg_type = Symbol
        default = :z
        help = "spin operators in the correlation function <ol or>"
    "--Ntau"
        arg_type = Int64
        default = 1601
        help = "number of τ"
    "--Niwn"
        arg_type = Int64
        default = 40
        help = "number of Masubara Frequency ωn"       
    "--bondD"
        arg_type = Int64
        default = 8
        help = "cMPS bond dimension"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/CMPO/TFIsing"
        help = "result folder"
end

parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const J = parsed_args[:J]
const Γ = parsed_args[:gamma]
const β = parsed_args[:beta]
const Ntau = parsed_args[:Ntau]
const Niwn = parsed_args[:Niwn]
const bondD = parsed_args[:bondD]

ResultFolder = parsed_args[:ResultFolder]
isdir(ResultFolder) || mkdir(ResultFolder)
ResultFolder = @sprintf "%s/J_%.2f_G_%.2f" ResultFolder J Γ
isdir(ResultFolder) || mkdir(ResultFolder)

const ol = parsed_args[:ol]
const or = parsed_args[:or]
#make_operators
opl = make_operator(pauli(ol), bondD)
opr = make_operator(pauli(or), bondD)

#CMPO
model = TFIsing(J, Γ)

CMPSPath = @sprintf "%s/CMPS/bondD_%02i_beta_%.2f.hdf5" ResultFolder bondD beta
ψ = readCMPS(CMPSPath)

tau = [i for i in range(0,β,length=Ntau)]
Gtau = [correlation_2time(t, opl, opr, ψ, model.Tmatrix, β) for t in tau]

GtauFolder= "$(ResultFolder)/Correlations/Gtau"
isdir(GtauFolder) || mkdir(GtauFolder)

GtauFile = @sprintf "%s/p%s_p%s_bondD_%02i_beta_%.2f.txt" GtauFolder ol or bondD β
open(GtauFile, "w") do file
    for i = 1:Ntau
        writedlm(file, [tau[i] Gtau[i]])
    end
end

wn = [Masubara_freq(n, β) for n=1:Niwn]
Giwn = [Masubara_freq_GF(n, opl, opr,ψ, model.Tmatrix, β) for n = 1:Niwn]

GiwnFolder= "$(ResultFolder)/Correlations/Giwn"
isdir(GiwnFolder) || mkdir(GiwnFolder)

GiwnFile = @sprintf "%s/p%s_p%s_bondD_%02i_beta_%.2f.txt" GiwnFolder ol or bondD β
open(GiwnFile, "w") do file
    for i = 1:Niwn
        writedlm(file, [wn[i] Giwn[i]])
    end
end


const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")

