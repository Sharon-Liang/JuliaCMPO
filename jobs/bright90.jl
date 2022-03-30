using Printf, DelimitedFiles
using Random; Random.seed!()

"""
    submitJob: Prepare a jobfile
"""
function submitJob(env, prog, args, jobname;
                  Run=false, 
                  Remove=true,
                  Wait=nothing, 
                  return_id = false,
                  proglang="julia")
    log_file_path = jobname * ".log"
    job =@sprintf """
    #!/bin/bash -l 
    #SBATCH --partition=a100 
    #SBATCH --nodes=1 
    #SBATCH --time=100:00:00 
    #SBATCH --cpus-per-task=4
    #SBATCH --job-name=%s 
    #SBATCH --output=%s 
    #SBATCH --error=%s 
    """ jobname log_file_path log_file_path
    
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
        sleep(2.0)
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