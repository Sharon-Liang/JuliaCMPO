using Printf
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

phys_model = "TFIsing"
env = "/home/sliang/JuliaCode/JuliaCMPO"
prog = env * "/jobs/CMPO_$(phys_model).jl"

bi = 11.3
bf = 40.0
bstep = 0.1
init = bi - bstep
Continue = true

Jlist = [1.0]
Γlist = [1.0]
init = 11.2
bondDlist = [16]



#CREAT LOG FOLDER
logdir = "/data/sliang/log/JuliaCMPO"
isdir(logdir) || mkdir(logdir)

for J in Jlist, Γ in Γlist, bondD in bondDlist
        args = Dict("J"=>J,
                    "G"=>Γ,
                    "bondD"=>bondD,
                    "bi"=>bi,
                    "bf"=>bf,
                    "bstep"=>bstep,
                    "init"=>init,
                    "Continue"=>Continue)
    jobname = logdir * "/" * phys_model
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/J_%.2f_G_%.2f_wid_01" jobname J Γ
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/bondD_%02i" jobname bondD
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/bi_%.2f_bf_%.2f_bstep_%.2f" jobname bi bf bstep

    jobid = submitJob(env, prog, args, jobname, 
                            Run = true)
end
