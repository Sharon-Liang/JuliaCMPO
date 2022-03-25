using Printf

function writeln(io::IO, xs...)
    write(io, xs...)
    write(io,"\n")
end

lambdaRange = [10^i for i in range(-6,0,length=13)]
betaRange = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0]
gRange =[1.0]
spectrum = [:A]
wmin = 1.e-8

#CREAT LOG FOLDER
logdir = "/home/sliang/JuliaCode/mycMPO/log/ising/nnls_dG"
isdir(logdir) || mkdir(logdir)

#CLEAR LOG FOLDER
#if length(readdir(logdir))!=0
#    for file in readdir(logdir)
#    run(```rm $(logdir)/$(file)```) end
#end

for s = 1:length(spectrum), j = 1:length(gRange)
    spec = @sprintf "%s" spectrum[s]
    g = @sprintf "%.1f" gRange[j]
    for b in 1:length(betaRange), i in 1:length(lambdaRange)
        beta = @sprintf "%.1f" betaRange[b]
        lambda = @sprintf "%e" lambdaRange[i]
        R = rand(Int)
        io = open("tmp$(R).sh","w+")
        writeln(io,"#!/bin/bash -l")
        writeln(io,"#SBATCH --partition=a100")
        writeln(io,"#SBATCH --nodes=1")
        writeln(io,"#SBATCH --time=999")
        writeln(io,"#SBATCH --job-name=NNLS_TFIsing_dG")
        writeln(io,"#SBATCH --output=$(logdir)/g_$(g)_beta_$(beta).log")
        writeln(io,"#SBATCH --error=$(logdir)/g_$(g)_beta_$(beta).log")
        writeln(io,"julia --project=/home/sliang/JuliaCode/mycMPO /home/sliang/JuliaCode/mycMPO/NNLS_TFIsing_dG.jl --g $(g) --beta $(beta) --lambda $(lambda) --spectrum $(spec) --wmin $(wmin)")
        close(io)
        println("Run: tmp$(R).sh")
        run(```sbatch tmp$(R).sh```)
        sleep(0.1)
        rm("tmp$(R).sh")
    end
end

    
