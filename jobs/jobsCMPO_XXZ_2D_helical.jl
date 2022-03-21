using Printf

beta_Range = [1.0, 16.0]
Jxy_Range = [1.0]
Jz_Range = [1.0]
bondD_Range = [10]
wid_Range = [1]

phys_model = "XXZ_2D_helical"

#CREAT LOG FOLDER
logdir = "/data/sliang/log/JuliaCMPO"
isdir(logdir) || mkdir(logdir)
logdir = "$(logdir)/$(phys_model)"
isdir(logdir) || mkdir(logdir)



for jxy = 1:length(Jxy_Range), jz = 1:length(Jz_Range), 
    d = 1:length(bondD_Range), w = 1:length(wid_Range),
    b = 1:length(beta_Range)

    Jxy = @sprintf "%.2f" Jxy_Range[jxy]
    Jz = @sprintf "%.2f" Jz_Range[jz]
    bondD = @sprintf "%02i" bondD_Range[d]
    width = @sprintf "%02i" wid_Range[w]
    β = @sprintf "%.2f" beta_Range[b]

    job_name = "Jz_$(Jz)_Jxy_$(Jxy)_wid_$(width)_D_$(bondD)_beta_$(β)" 
    log_file_path = "$(logdir)/$(job_name).log"
    
    R = rand(Int)
    io = open("tmp$(R).sh","w+")
    write(io,"#!/bin/bash -l \n\
            #SBATCH --partition=a100 \n\
            #SBATCH --cpus-per-task=1 \n\
            #SBATCH --time=999 \n\
            #SBATCH --job-name=$(job_name) \n\
            #SBATCH --output=$(log_file_path) \n\
            #SBATCH --error=$(log_file_path) \n\
            julia --project=/home/sliang/JuliaCode/JuliaCMPO \
                /home/sliang/JuliaCode/JuliaCMPO/jobs/CMPO_$(phys_model).jl \
                --Jxy $(Jxy) --Jz $(Jz) --width $(width) \
                --beta $(β) --bondD $(bondD)"
            )
    close(io)
    println("Run: tmp$(R).sh")
    run(```sbatch tmp$(R).sh```)
    sleep(0.1)
    rm("tmp$(R).sh")
end

# \n\ 换行顶格，后一个\之后不能有任何空格