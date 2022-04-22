using Printf, Dates
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

phys_model = "XXZ_2D_helical"
env = "/home/sliang/JuliaCode/JuliaCMPO"
prog = env * "/jobs/CMPO_$(phys_model).jl"

device = "gpu"
logtag = Dates.format(now(), "yyyy-mm-dd")*"-"*device

Wait = nothing
cpu_per_task = 8

tag = "2022-04-20"*"-"*device
beta1 = [0.1, 0.5]
beta2 = [i for i in range(1.0, 10.0, step=1.0)]
beta3 = [i for i in range(1.1, 1.9, step=0.1)]

#βlist = [0.1, 0.5]
βlist = vcat(beta1, beta2, beta3)
Jzlist = [1.0]
Jxylist = [1.0]
bondDlist = [32]
widlist = [3, 5]
Continue = 999 #Continue > max_pow_step,  Continue = true
max_pow_step = 100
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
                            machine = "a100",
                            usegpu = usegpu,
                            cpu_per_task = cpu_per_task,
                            Run = true, 
                            Wait = Wait)
    end   
end      
