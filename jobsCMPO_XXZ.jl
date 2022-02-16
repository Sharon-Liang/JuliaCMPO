using Printf

function writeln(io::IO, xs...)
    write(io, xs...)
    write(io,"\n")
end

beta_max = 30.0
#beta_min = 1.0
#beta_step = 0.1  #defalt β-format %.2f
JzRange =[1.0]
WidthRange = [3.0]


#CREAT LOG FOLDER
logdir = "/home/sliang/JuliaCode/mycMPO/log/xxz"
isdir(logdir) || mkdir(logdir)
logdir = logdir*"/cmpo"
isdir(logdir) || mkdir(logdir)

#CLEAR LOG FOLDER
#if length(readdir(logdir))!=0
#    for file in readdir(logdir)
#    run(```rm $(logdir)/$(file)```) end
#end

for j = 1:length(JzRange), w = 1:length(WidthRange)
    Jz = @sprintf "%.1f" JzRange[j]
    width = @sprintf "%i" WidthRange[w]
    R = rand(Int)
    io = open("tmp$(R).sh","w+")
    writeln(io,"#!/bin/bash -l")
    writeln(io,"#SBATCH --partition=a100")
    writeln(io,"#SBATCH --nodes=1")
    writeln(io,"#SBATCH --time=999")
    writeln(io,"#SBATCH --job-name=CMPO_XXZ")
    writeln(io,"#SBATCH --output=$(logdir)/Jz_$(Jz)_W_$(width).log")
    writeln(io,"#SBATCH --error=$(logdir)/Jz_$(Jz)_W_$(width).log")
    writeln(io,"julia --project=/home/sliang/JuliaCode/mycMPO /home/sliang/JuliaCode/mycMPO/CMPO_XXZ.jl --Jz $(Jz) --width $(width) --beta_max $(beta_max)")
    close(io)
    println("Run: tmp$(R).sh")
    run(```sbatch tmp$(R).sh```)
    sleep(0.1)
    rm("tmp$(R).sh")
end


    
