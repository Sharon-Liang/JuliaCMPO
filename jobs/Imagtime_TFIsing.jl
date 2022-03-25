#using Pkg; Pkg.activate("/home/sliang/JuliaCode/mycMPO/")
using cMPO
using TimerOutputs, Dates, ArgParse
using DelimitedFiles, JLD, HDF5, Printf
using LinearAlgebra

const to = TimerOutput()
const Start_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

settings = ArgParseSettings(prog="G(τ) and dGtau for TFIsing model"
)
@add_arg_table! settings begin
    "--g"
        arg_type = Float64
        default = 1.0
        help = "Transverse Field Ising Γ/J"
    "--beta"
        arg_type = Float64
        default = 10.0
        help = "β: inverse temperature"
    "--D"
        arg_type = Int64
        default = 8
        help = "cMPO bond dimension"
    "--operator"
        arg_type = Symbol
        default = :z
        help = "spin operators in the correlation function"
    "--Ntau"
        arg_type = Int64
        default = 1601
        help = "number of τ"
    "--DataFolder"
        arg_type = String
        default = "/data/sliang/CMPO/ising"
        help = "data folder"
    "--ResultFolder"
        arg_type = String
        default = "/data/sliang/CMPO/ising/imagtime"
        help = "result folder"
    "--dodGtau"
        arg_type = Bool
        default = true
        help = "calculate dG(τ) or not"
end
parsed_args = parse_args(settings; as_symbols=true)
print(parsed_args,"\n")

const g = parsed_args[:g]
const β = parsed_args[:beta]
const D = parsed_args[:D]
const operator = parsed_args[:operator]
const Ntau = parsed_args[:Ntau]
const DataFolder = parsed_args[:DataFolder]
const ResultFolder = parsed_args[:ResultFolder]
const dodGtau = parsed_args[:dodGtau]

#DataFile Path 
DataFile1 = @sprintf "%s/D_%i/g_%.1f.jld" DataFolder D g

#Existence of DataFile
@assert isfile(DataFile1)

#ResultFile Path
Gtau_ResultFile1 = @sprintf "%s/gtau/g_%.1f_D_%i_beta_%i.txt" ResultFolder g D β
Gtau_ResultFile2 = @sprintf "%s/gtau/g_%.1f_D_%im2_beta_%i.txt" ResultFolder g D β
#Remove old ResultFile in the first step
isfile(Gtau_ResultFile1) && rm(Gtau_ResultFile1) 

#LOAD DATA
@timeit to "load data" begin
    d1 = load(DataFile1)
end
w = TFIsing(1.0, g)
o = pauli(operator)
o1 = make_operator(o, D)
o2 = make_operator(o, 2D)

ψ1 = tocmps(d1[string(β)][2])
ψ2 = w * ψ1

tau = [i for i in range(0,β,length=Ntau)]
@timeit to "G(τ)" begin
    Gt1 = [correlation_2time(τ,o1,o1',ψ1,w,β) for τ in tau]
    Gt2 = [correlation_2time(τ,o2,o2',ψ2,w,β) for τ in tau]
end

#WRITE RESULT
open(Gtau_ResultFile1,"w+") do file
    for i=1: Ntau
        writedlm(file,[tau[i] Gt1[i]])
    end
end

open(Gtau_ResultFile2,"w+") do file
    for i=1: Ntau
        writedlm(file,[tau[i] Gt2[i]])
    end
end


if dodGtau
    if operator == :z
        õ = 2g .* pauli(:iy)
    else
        @error("not support yet")
    end

    dGtau_ResultFile1 = @sprintf "%s/dgtau/g_%.1f_D_%i_beta_%i.txt" ResultFolder g D β
    dGtau_ResultFile2 = @sprintf "%s/dgtau/g_%.1f_D_%im2_beta_%i.txt" ResultFolder g D β
    isfile(dGtau_ResultFile1) && rm(dGtau_ResultFile1) 

    õ1 = make_operator(õ, D)
    õ2 = make_operator(õ, 2D)

    @timeit to "dG(τ)" begin
        dGt1 = [correlation_2time(τ,õ1,o1',ψ1,w,β) for τ in tau]
        dGt2 = [correlation_2time(τ,õ2,o2',ψ2,w,β) for τ in tau]
    end

    open(dGtau_ResultFile1,"w+") do file
        for i=1: Ntau
            writedlm(file,[tau[i] dGt1[i]])
        end
    end
    
    open(dGtau_ResultFile2,"w+") do file
        for i=1: Ntau
            writedlm(file,[tau[i] dGt2[i]])
        end
    end
end

const End_Time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
const Running_TimeTable = string(to)
@show Start_Time
@show End_Time
print(Running_TimeTable,"\n")



