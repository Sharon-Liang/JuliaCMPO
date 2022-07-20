using Printf, Dates
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

phys_model = "XXZ_2D_helical"
env = "/home/sliang/JuliaCode/JuliaCMPO"
prog = env * "/jobs/CMPO_$(phys_model).jl"

processor = GPU
machine = p100
gpu_memory = 80
logtag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"

Wait = nothing
cpu_per_task = 1

#tag = "2022-06-21"*"-"*processor
tag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"

βlist = [0.5] 
Jzlist = [1.0]
Jxylist = [1.0]
bondDlist = [32]
widlist = [3]
Continue = 0  #Continue > max_pow_step,  Continue = true
max_pow_step = 10



#CREAT LOG FOLDER
logdir = "/data/sliang/JuliaCMPO/log"
isdir(logdir) || mkdir(logdir)

for bondD in bondDlist
    for Jxy in Jxylist, Jz in Jzlist, width in widlist, β in βlist
        args = Dict("Jz"=>Jz,
                    "Jxy"=>Jxy,
                    "bondD"=>bondD,
                    "width"=>width,
                    "beta"=>β,
                    "max_pow_step"=>max_pow_step,
                    "Continue"=>Continue,
                    "tag"=>tag,
                    "processor"=> Int(processor)
                    )
        jobname = logdir * "/" * phys_model
        isdir(jobname) || mkdir(jobname)
        jobname = @sprintf "%s/Jz_%.2f_Jxy_%.2f_wid_%02i" jobname Jz Jxy width
        isdir(jobname) || mkdir(jobname)
        jobname = @sprintf "%s/bondD_%02i_%s" jobname bondD logtag
        isdir(jobname) || mkdir(jobname)
        timetag = Dates.format(now(), "HH-MM-SS")
        jobname = @sprintf "%s/beta_%.2f_%s" jobname β timetag

        jobid = submitJob(env, prog, args, jobname, 
                            partitian = machine,
                            processor = processor,
                            gpu_memory = gpu_memory,
                            cpu_per_task = cpu_per_task,
                            Run = true, 
                            Wait = Wait)
    end   
end      
