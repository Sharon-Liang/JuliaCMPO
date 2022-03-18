using TimerOutputs, Dates, ArgParse
using DelimitedFiles, JLD, HDF5, Printf
using cMPO
using Optim

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
settings = ArgParseSettings(prog="Correlations: TFIsing model"
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
    "--beta"
        arg_type = Float64
        default = 1.0
        help = "β: inverse temperature"
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
        default = "/data/sliang/CMPO/XXZ_2D_helical"
        help = "result folder"
end

parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const Jxy = parsed_args[:Jxy]
const Jz = parsed_args[:Jz]
const width = parsed_args[:width]
const β = parsed_args[:beta]
const Ntau = parsed_args[:Ntau]
const Niwn = parsed_args[:Niwn]
const bondD = parsed_args[:bondD]


ResultFolder = parsed_args[:ResultFolder]
isdir(ResultFolder) || mkdir(ResultFolder)
ResultFolder = @sprintf "%s/Jxy_%.2f_Jz_%.2f_wid_%02i" ResultFolder Jxy Jz width
isdir(ResultFolder) || mkdir(ResultFolder)

model = XXZmodel_2D_helical(Jz/Jxy, width)
CMPSPath = @sprintf "%s/CMPS/bondD_%02i_beta_%.2f.hdf5" ResultFolder bondD β
ψ = readCMPS(CMPSPath)


operator = :z
op_name = "sz_sz"
op = make_operator(0.5 * pauli(:z), 8)


ResultFolder = "$(ResultFolder)/Correlations"
isdir(ResultFolder) || mkdir(ResultFolder)

tau = [i for i in range(0,β,length=Ntau)]
@timeit to "Gtau" begin
    Gtau = [correlation_2time(t, op, op, ψ, model.Tmatrix, β) for t in tau]
end

GtauFolder= "$(ResultFolder)/Gtau"
isdir(GtauFolder) || mkdir(GtauFolder)

GtauFile = @sprintf "%s/%s_bondD_%02i_beta_%.2f.txt" GtauFolder op_name bondD β
open(GtauFile, "w") do file
    for i = 1:Ntau
        writedlm(file, [tau[i] Gtau[i]])
    end
end

wn = [Masubara_freq(n, β) for n=1:Niwn]
@timeit to "Giωn" begin
    Giwn = [Masubara_freq_GF(n, op, op, ψ, model.Tmatrix, β) for n = 1:Niwn]
end

GiwnFolder= "$(ResultFolder)/Giwn"
isdir(GiwnFolder) || mkdir(GiwnFolder)

GiwnFile = @sprintf "%s/%s_bondD_%02i_beta_%.2f.txt" GiwnFolder op_name bondD β
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

