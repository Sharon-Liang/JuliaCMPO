using Printf
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

phys_model = "TFIsing_2D_helical"
env = "/home/sliang/JuliaCode/JuliaCMPO"
prog = env * "/jobs/CMPO_$(phys_model).jl"

processor = GPU
machine = a100
gpu_memory = 80
logtag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"

Wait = nothing
cpu_per_task = 8

#tag = "2022-06-21"*"-"*processor
tag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"

βlist = [0.1]
Jlist = [1.0]
Γlist = [2.0]
bondDlist = [8]
widlist = [4]
Continue = 0  #Continue > max_pow_step,  Continue = true
max_pow_step = 100

#CREAT LOG FOLDER
logdir = "/data/sliang/log/JuliaCMPO"
isdir(logdir) || mkdir(logdir)

for bondD in bondDlist
    for J in Jlist, Γ in Γlist, width in widlist, β in βlist
        args = Dict("J"=>J,
                    "G"=>Γ,
                    "bondD"=>bondD,
                    "width"=>width,
                    "beta"=>β,
                    "max_pow_step"=>max_pow_step,
                    "Continue"=>Continue,
                    "tag"=>tag,
                    "processor"=>processor
                    )
        jobname = logdir * "/" * phys_model
        isdir(jobname) || mkdir(jobname)
        jobname = @sprintf "%s/J_%.2f_G_%.2f_wid_%02i_unexpanded" jobname J Γ width
        isdir(jobname) || mkdir(jobname)
        jobname = @sprintf "%s/bondD_%02i" jobname bondD
        isdir(jobname) || mkdir(jobname)
        jobname = @sprintf "%s/beta_%.2f" jobname β

        jobid = submitJob(env, prog, args, jobname, 
                            partitian = machine,
                            processor = processor,
                            gpu_memory = gpu_memory,
                            cpu_per_task = cpu_per_task,
                            Run = true, 
                            Wait = Wait)
    end   
end      
