#module evaluate

"""
    evaluate(A, bondD, β, ResultFolder; evaluate_options)
    
Calculate left and right CMPS of a CMPO `A` at temperature `β`,
results are save in folder `ResultFolder` 
"""
function evaluate(A::AbstractCMPO, bondD::Integer, β::Real, ResultFolder::String; evaluate_options::EvaluateOptions)
    @unpack hermitian = evaluate_options
    hermitian === nothing ? hermitian = ishermitian(A) : hermitian = hermitian
    if hermitian
        hermitian_evaluate(A, bondD, β, ResultFolder; evaluate_options)
    else
        non_hermitian_evaluate(A, bondD, β, ResultFolder; evaluate_options)
    end
end


"""
    hermitian_evaluate(A, bondD, β, ResultFolder; evaluate_options)

Calculate left and right CMPS of a CMPO `A` at temperature `β` by optimizing free_energy,
results are save in folder `ResultFolder` 
"""
function hermitian_evaluate(A::AbstractCMPO, bondD::Integer, β::Real, ResultFolder::String; 
                            evaluate_options::EvaluateOptions)
    @unpack (init, processor, trace_estimator, 
             compress_options, optim_options, tag, show_trace) = evaluate_options
    solver = solver_function(processor)

    """
    ResultFolder: classified by Model parameters(interaction, width), store CMPS and Obsv files
    CMPSResultFolder: CMPS information, classified by bond dimension
    """
    trace_estimator === nothing ? estimator_name = "nothing" : estimator_name = string(trace_estimator.estimator)
    CMPSResultFolder = @sprintf "%s/bondD_%02i_CMPS_%s_%s" ResultFolder bondD estimator_name tag
    OptResultFolder = @sprintf "%s/bondD_%02i_Opt_%s_%s" ResultFolder bondD estimator_name tag
    isdir(ResultFolder) || mkdir(ResultFolder)
    isdir(CMPSResultFolder) || mkdir(CMPSResultFolder) 
    isdir(OptResultFolder) || mkdir(OptResultFolder)

    A = solver(A)
    #initiate cmps
    init === nothing ? 
        ψ = init_cmps(bondD, A, β; compress_options) : ψ = solver(init)

    ψ = diagQ(ψ)

    dQ = convert(Vector, diag(ψ.Q))
    R  = convert(Array, ψ.R)

    loss() = free_energy(solver(cmps_generate(diagm(dQ), R)), A, β, trace_estimator)
    F_initial = loss()

    p0, f, g! = optim_functions(loss, Zygote.Params([dQ, R]))
    optim_result = Optim.optimize(f, g!, p0, LBFGS(), optim_options)
    F_final = minimum(optim_result)

    ψ = cmps_generate(diagm(dQ), R) |> solver
    #save optimize result 
    OptResultFile = @sprintf "%s/beta_%.2f.txt" OptResultFolder β
    open(OptResultFile, "w") do file
        F_initial_string = @sprintf "F_initial = %.16f \n" F_initial
        F_final_string =   @sprintf "F_final   = %.16f \n" F_final
        write(file, F_initial_string)
        write(file, F_final_string)
        write(file, optim_result)
        write(file, optim_result.trace)
    end

    # calculate thermal dynamic quanties
    dict = Dict()
    dict["F"] = F_final
    #dict["E"] = energy(ψ, A, β, trace_estimator)
    #dict["Cv"] = specific_heat(ψ, A, β)
    #dict["S"] = β * (dict["E"] - dict["F"])
    ResultFile = @sprintf "%s/beta_%.2f.hdf5" CMPSResultFolder β
    saveCMPS(ResultFile, CTensor(ψ), dict)
    return  ψ, dict
end


"""
    non_hermitian_evaluate(A, bondD, β, ResultFolder; evaluate_options)

Calculate left and right CMPS of a CMPO `A` at temperature `β` by power method,
results are save in folder `ResultFolder` 
"""
function non_hermitian_evaluate(A::AbstractCMPO, bondD::Integer, β::Real, ResultFolder::String; 
                                evaluate_options::EvaluateOptions)
    @unpack (init, processor, trace_estimator, compress_options, 
             tocontinue, max_pow_step, group, tag, show_trace) = evaluate_options
    solver = solver_function(processor)
    
    """
        ResultFolder: classified by Model parameters(interaction, width), store CMPS and Obsv files
        CMPSResultFolder: CMPS information, classified by bond dimension
        ChkpFolder: Check points
    """
    A = solver(A)
    g = 1
    while g < group #enlarged matrix: phy_dim -> phy_dim^group
        A *= A
        g += 1
    end

    trace_estimator === nothing ? estimator_name = "nothing" : estimator_name = string(trace_estimator.estimator)
    CMPSResultFolder = @sprintf "%s/bondD_%02i_CMPS_%s_%s" ResultFolder bondD estimator_name tag
    ChkpFolder = @sprintf "%s/CheckPoint_beta_%.2f" CMPSResultFolder β
    ChkpPsiFolder = @sprintf "%s/cmps_psi" ChkpFolder
    ChkpLpsiFolder = @sprintf "%s/cmps_Lpsi" ChkpFolder
    ChkpOptPsiFolder = @sprintf "%s/opt_psi" ChkpFolder
    ChkpOptLpsiFolder = @sprintf "%s/opt_Lpsi" ChkpFolder
    isdir(ResultFolder) || mkdir(ResultFolder)
    isdir(CMPSResultFolder) || mkdir(CMPSResultFolder) 
    isdir(ChkpFolder) || mkdir(ChkpFolder)

    ChkpEngFile = "$(ChkpFolder)/Obsv_FECvS.txt"
    PsiFidelityFile = "$(ChkpFolder)/Fidelity_psi.txt"
    LpsiFidelityFile = "$(ChkpFolder)/Fidelity_Lpsi.txt"

    if tocontinue == false || tocontinue == 0
        #CLEAR ChkpFolder
        if length(readdir(ChkpFolder)) != 0
            for file in readdir(ChkpFolder) rm("$(ChkpFolder)/$(file)",recursive=true) end
        end
        isdir(ChkpPsiFolder) || mkdir(ChkpPsiFolder)
        isdir(ChkpLpsiFolder) || mkdir(ChkpLpsiFolder)
        isdir(ChkpOptPsiFolder) || mkdir(ChkpOptPsiFolder)
        isdir(ChkpOptLpsiFolder) || mkdir(ChkpOptLpsiFolder)

        #initiate cmps
        pow_step = 0
        if init === nothing 
            ψr = init_cmps(bondD, A, β; compress_options)
            ψl = init_cmps(bondD, transpose(A), β; compress_options)
        else 
            ψr, ψl = solver(init)
        end

        ψr, ψl = diagQ(ψr), diagQ(ψl) 
        ChkpPsiFile = @sprintf "%s/step_%03i.hdf5" ChkpPsiFolder pow_step
        ChkpLpsiFile = @sprintf "%s/step_%03i.hdf5" ChkpLpsiFolder pow_step
        saveCMPS(ChkpPsiFile, CTensor(ψr))
        saveCMPS(ChkpLpsiFile, CTensor(ψl))

        open(ChkpEngFile,"w") do cfile
            F = free_energy(ψl, ψr, A, β)/group
            #E = energy(ψl, ψr, A, β)/group
            #Cv = specific_heat(ψl, ψr, A, β)/group
            #S = β * (E - F)    
            #write(cfile, "step    free_energy/site        energy/site         entropy/site    \n")
            #write(cfile, "----  -------------------  --------------------  -------------------\n")
            #EngString = @sprintf "%3i   %.16f   %.16f   %.16f \n" pow_step F E S
            write(cfile, "step    free_energy/site   \n")
            write(cfile, "----  -------------------  \n")
            EngString = @sprintf "%3i   %.16f \n" pow_step F 
            write(cfile, EngString)
        end
        for pfile in [PsiFidelityFile, LpsiFidelityFile]
            open(pfile,"w") do cfile
                write(cfile, "step   initial fidelity      final fidelity  \n")
                write(cfile, "----  ------------------   ------------------\n")
            end
        end
    else
        if tocontinue == true
            pow_step = readdlm(ChkpEngFile, skipstart=2)[end,1]
        else
            pow_step = tocontinue
            EngData = readlines(ChkpEngFile, keep=true)
            open(ChkpEngFile, "w") do file
                for i = 1:pow_step+3 write(file, EngData[i]) end
            end
            for pfile in [PsiFidelityFile, LpsiFidelityFile]
                psi_fidelity = readlines(pfile, keep=true)
                open(pfile, "w") do file
                    for i = 1:pow_step+2 write(file, psi_fidelity[i]) end
                end
            end
        end
        ChkpPsiFile = @sprintf "%s/step_%03i.hdf5" ChkpPsiFolder pow_step
        ChkpLpsiFile = @sprintf "%s/step_%03i.hdf5" ChkpLpsiFolder pow_step
        ψr = readCMPS(ChkpPsiFile)
        ψl = readCMPS(ChkpLpsiFile)
    end

    while pow_step < max_pow_step
        pow_step += 1
        if show_trace println(@sprintf "Power Step = %03i" pow_step) end

        ψr, ψl = solver(ψr), solver(ψl)
        Psi = compress_cmps(A * ψr, bondD, β; compress_options)
        Lpsi = compress_cmps(transpose(A) * ψl, bondD, β; compress_options)
        ψr, ψl = Psi.cmps, Lpsi.cmps

        open(ChkpEngFile,"a") do cfile
            F = free_energy(ψl, ψr, A, β)/group
            #E = energy(ψl, ψr, A, β)/group
            #Cv = specific_heat(ψl, ψr, A, β)/group
            #S = β * (E - F)
            #EngString = @sprintf "%3i   %.16f   %.16f   %.16f \n" pow_step F E S
            EngString = @sprintf "%3i   %.16f \n" pow_step F 
            write(cfile, EngString)
        end

        ChkpPsiFile = @sprintf "%s/step_%03i.hdf5" ChkpPsiFolder pow_step
        ChkpLpsiFile = @sprintf "%s/step_%03i.hdf5" ChkpLpsiFolder pow_step
        saveCMPS(ChkpPsiFile, CTensor(ψr))
        saveCMPS(ChkpLpsiFile, CTensor(ψl))

        open(PsiFidelityFile,"a") do cfile
            FidelityString = @sprintf "%3i   %.16f   %.16f \n" pow_step Psi.ifidel Psi.ffidel
            write(cfile, FidelityString)
        end

        OptPsitFile = @sprintf "%s/step_%03i.txt" ChkpOptPsiFolder pow_step
        open(OptPsitFile, "w") do file
            write(file, Psi.optim_result)
            if Psi.optim_result.trace !== nothing write(file, Psi.optim_result.trace) end
        end

        open(LpsiFidelityFile,"a") do cfile
            FidelityString = @sprintf "%3i   %.16f   %.16f \n" pow_step Lpsi.ifidel Lpsi.ffidel
            write(cfile, FidelityString)
        end

        OptLpsiFile = @sprintf "%s/step_%03i.txt" ChkpOptLpsiFolder pow_step
        open(OptLpsiFile, "w") do file
            write(file, Lpsi.optim_result)
            if Lpsi.optim_result.trace !== nothing write(file, Lpsi.optim_result.trace) end
        end
    end
    return 
end

#end    # module evaluate