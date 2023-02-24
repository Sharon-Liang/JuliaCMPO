using Printf, Dates
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

phys_model = "J1J2Chain"
env = "/home/sliang/JuliaCode/JuliaCMPO"
prog = env * "/jobs/CMPO_$(phys_model).jl"

processor = GPU
machine = titanv

logtag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"
Wait = nothing
cpu_per_task = 8

tag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"

J1list = [1.0]
J2list = [0.241167, 0.5]
bondDlist = [20]
βlist = [20, 30]
to_group = 2
to_shift = 0.


#CREAT LOG FOLDER
logdir = "/data/sliang/JuliaCMPO/log"
isdir(logdir) || mkdir(logdir)

for J1 in J1list, J2 in J2list, bondD in bondDlist, β in βlist
        args = Dict("J1"=>J1,
                    "J2"=>J2,
                    "bondD"=>bondD,
                    "beta"=>β,
                    "tag"=>tag,
                    "group"=>to_group,
                    "shift"=>to_shift,
                    "processor"=> Int(processor)
                    )
    jobname = logdir * "/" * phys_model
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/J1_%f_J2_%f_wid_01_bondD_%02i_group_%i_shift_%e_%s" jobname J1 J2 bondD to_group to_shift logtag
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/beta_%f" jobname β

    jobid = submitJob(env, prog, args, jobname, 
                            partitian = machine,
                            processor = processor,
                            cpu_per_task = cpu_per_task,
                            Run = false,
                            Remove=false, 
                            Wait = Wait)
end
