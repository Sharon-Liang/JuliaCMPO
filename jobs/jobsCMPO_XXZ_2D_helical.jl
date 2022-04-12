using Printf, Dates
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

phys_model = "XXZ_2D_helical"
env = "/home/sliang/JuliaCode/JuliaCMPO"
prog = env * "/jobs/CMPO_$(phys_model).jl"

device = "gpu"
logtag = Dates.format(now(), "yyyy-mm-dd")*"-"*device

Wait = nothing
cpu_per_task = 4

tag = logtag
βlist = [1.0]
Jzlist = [1.0]
Jxylist = [1.0]
bondDlist = [16]
widlist = [5]
Continue = 0 #Continue > max_pow_step,  Continue = true
max_pow_step = 2
device == "gpu" ? usegpu = true : usegpu = false



#CREAT LOG FOLDER
logdir = "/data/sliang/log/JuliaCMPO"
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
                    "device"=> device
                    )
        jobname = logdir * "/" * phys_model
        isdir(jobname) || mkdir(jobname)
        jobname = @sprintf "%s/Jz_%.2f_Jxy_%.2f_wid_%02i" jobname Jz Jxy width
        isdir(jobname) || mkdir(jobname)
        jobname = @sprintf "%s/bondD_%02i_%s" jobname bondD logtag
        isdir(jobname) || mkdir(jobname)
        jobname = @sprintf "%s/beta_%.2f" jobname β

        jobid = submitJob(env, prog, args, jobname, 
                            machine = "p100",
                            usegpu = usegpu,
                            cpu_per_task = cpu_per_task,
                            Run = true, 
                            Wait = Wait)
    end   
end      
