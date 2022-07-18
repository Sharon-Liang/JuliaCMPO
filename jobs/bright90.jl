using Printf, DelimitedFiles
using Random; Random.seed!()

@enum Partitian a100 p100 v100 titanv
@enum Processor CPU GPU

"""
    submitJob: Prepare a jobfile
"""
function submitJob(env, prog, args, jobname;
                  partitian::Partitian = a100,
                  processor::Processor = CPU,
                  gpu_memory::Int = 0,
                  cpu_per_task::Integer = 4,
                  Run=false, 
                  Remove=true,
                  Wait=nothing, 
                  return_id = false,
                  proglang="julia",
                  runtime::String = "99-99:99:99")
    log_file_path = jobname * ".log"
    job =@sprintf """
    #!/bin/bash -l 
    #SBATCH --partition=%s 
    #SBATCH --nodes=1 
    #SBATCH --time=%s 
    #SBATCH --cpus-per-task=%i
    #SBATCH --job-name=%s 
    #SBATCH --output=%s 
    #SBATCH --error=%s 
    """ partitian runtime cpu_per_task jobname log_file_path log_file_path
    
    if device == GPU
        if partitian == a100 && gpu_memory != 0
            job *= """
            #SBATCH --gres=gpu:A100_$(gpu_memory)G:1
            """
        else 
            job *= """
            #SBATCH --gres=gpu:1
            """
        end
    end

    if Wait !== nothing
        dependency = @sprintf """
        #SBATCH --dependency=afterany:%d
        """ Wait
        job *= dependency
    end
    
    if return_id
        job *= """
        echo "\$SLURM_JOB_ID" > $(jobname)_id.txt
        """
    end
    
    job *="""
    echo "The current job ID is \$SLURM_JOB_ID"
    echo "Running on \$SLURM_JOB_NUM_NODES nodes:"
    echo \$SLURM_JOB_NODELIST
    echo "Using \$SLURM_CPUS_PER_TASK cpus-per-task"
    echo "Using \$SLURM_NTASKS_PER_NODE tasks per node"
    echo "A total of \$SLURM_NTASKS tasks is used"
    echo "CUDA devices \$CUDA_VISIBLE_DEVICES"
    """

    job *= """
    echo Job started at `date`
    """

    job *= "$(proglang) --threads \$SLURM_CPUS_PER_TASK --project=$(env) $(prog) "

    for (key, val) in args
        job *= @sprintf "--%s %s " key val
    end

    job *= """
    \necho Job finished at `date`
    """

    R = rand(Int)
    jobfile = "tmp$(R).sh"
    open(jobfile, "w") do file
        write(file, job)
    end

    jobid = nothing

    if Run 
        run(`sbatch $(jobfile)`) 
        sleep(1.0)
        if return_id
            jobid = readdlm(jobname*"_id.txt", Int64)[1]
            rm(jobname*"_id.txt")
        end
    else
        run(`cat $(jobfile)`)
    end
    if Remove rm(jobfile) end
    return jobid
end