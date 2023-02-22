using Printf, Dates
include("/home/sliang/JuliaCode/JuliaCMPO/jobs/bright90.jl")

phys_model = "J1J2Chain"
env = "/home/sliang/JuliaCode/JuliaCMPO"
prog = env * "/jobs/CMPO_$(phys_model).jl"

processor = CPU
machine = titanv

logtag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"
Wait = nothing
cpu_per_task = 8

tag = Dates.format(now(), "yyyy-mm-dd")*"-"*"$(processor)"

J1list = [1.0]
J2list = [0.5]
bondDlist = [20]


#CREAT LOG FOLDER
logdir = "/data/sliang/JuliaCMPO/log"
isdir(logdir) || mkdir(logdir)

for J1 in J1list, J2 in J2list, bondD in bondDlist
        args = Dict("J1"=>J1,
                    "J2"=>J2,
                    "bondD"=>bondD,
                    "tag"=>tag,
                    "processor"=> Int(processor)
                    )
    jobname = logdir * "/" * phys_model
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/J1_%f_J2_%f_wid_01_bondD_%02i_%s" jobname J1 J2 bondD logtag
    isdir(jobname) || mkdir(jobname)
    jobname = @sprintf "%s/bi_%f_bf_%f" jobname 1/0.05 1/0.03 

    jobid = submitJob(env, prog, args, jobname, 
                            partitian = machine,
                            processor = processor,
                            cpu_per_task = cpu_per_task,
                            Run = true, 
                            Wait = Wait)
end
