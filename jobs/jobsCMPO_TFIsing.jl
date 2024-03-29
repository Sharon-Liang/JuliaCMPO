using Printf, Dates
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

phys_model = "TFIsing"
env = "/home/sliang/JuliaCode/JuliaCMPO"
prog = env * "/jobs/CMPO_$(phys_model).jl"

processor = GPU
machine = titanv

logtag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"
Wait = nothing
cpu_per_task = 8

tag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"

bi = 1.0
bf = 1.5
bstep = 0.1

Jlist = [1.0]
Γlist = [1.0]
bondDlist = [8]


#CREAT LOG FOLDER
logdir = "/data/sliang/JuliaCMPO/log"
isdir(logdir) || mkdir(logdir)

for J in Jlist, Γ in Γlist, bondD in bondDlist
        args = Dict("J"=>J,
                    "G"=>Γ,
                    "bondD"=>bondD,
                    "bi"=>bi,
                    "bf"=>bf,
                    "bstep"=>bstep,
                    "tag"=>tag,
                    "processor"=> Int(processor)
                    )
    jobname = logdir * "/" * phys_model
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/J_%.2f_G_%.2f_wid_01_bondD_%02i_%s" jobname J Γ bondD logtag
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/bi_%.2f_bf_%.2f_bstep_%.2f" jobname bi bf bstep

    jobid = submitJob(env, prog, args, jobname, 
                            partitian = machine,
                            processor = processor,
                            cpu_per_task = cpu_per_task,
                            Run = true, 
                            Wait = Wait)
end
