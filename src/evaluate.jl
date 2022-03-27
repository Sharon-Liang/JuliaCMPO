#module evaluate
"""
    Evaluate PhysModel m when the hermiticty of its transfer matrix is unknown or not specified
"""
function evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
                    init::Union{CMPS, Nothing} = nothing, 
                    max_pow_step::Integer = 100,
                    hermitian::Union{Bool, Nothing} = nothing,
                    show_trace::Bool = false)
    hermitian === nothing ? hermitian = ishermitian(m.Tmatrix) : hermitian = hermitian
    if hermitian
        hermitian_evaluate(m, bondD, β, ResultFolder, 
            init = init, 
            show_trace = show_trace)
    else
        non_hermitian_evaluate(m, bondD, β, ResultFolder, 
            init = init, 
            max_pow_step = max_pow_step,
            show_trace = show_trace)
    end
end


"""
    Evaluate PhysModel m when its transfer matrix is hermitian
"""
function hermitian_evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
            init::Union{CMPS, Nothing} = nothing,
            show_trace::Bool=false)
    """
    ResultFolder: classified by Model parameters(interaction, width), store CMPS and Obsv files
    CMPSResultFolder: CMPS information, classified by bond dimension
    """
    CMPSResultFolder = @sprintf "%s/bondD_%02i_CMPS" ResultFolder bondD
    OptResultFolder = @sprintf "%s/bondD_%02i_Opt" ResultFolder bondD
    isdir(ResultFolder) || mkdir(ResultFolder)
    isdir(CMPSResultFolder) || mkdir(CMPSResultFolder) 
    isdir(OptResultFolder) || mkdir(OptResultFolder)

    Tm = m.Tmatrix
    #initiate cmps
    init === nothing ? 
        ψ = init_cmps(bondD, Tm, β, show_trace = show_trace) : ψ = init

    ψ = ψ |> diagQ
    loss() = free_energy(CMPS(diag(ψ.Q)|> diagm, ψ.R), Tm, β)
    F_initial = loss()

    p0, f, g! = optim_functions(loss, Params([ψ.Q, ψ.R]))
    optim_options = Optim.Options(f_tol = eps(), g_tol = 1.e-8,
                            iterations = 10000,
                            store_trace = true,
                            show_trace = show_trace, show_every = 10)
    optim_result = Optim.optimize(f, g!, p0, LBFGS(), optim_options)
    F_final = minimum(optim_result)

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
    dict["E"] = energy(ψ, Tm, β)
    dict["Cv"] = specific_heat(ψ, Tm, β)
    dict["S"] = β * (dict["E"] - dict["F"])
    ResultFile = @sprintf "%s/beta_%.2f.hdf5" CMPSResultFolder β
    saveCMPS(ResultFile, ψ, dict)
    return  ψ, dict
end


"""
    Evaluate PhysModel m when its transfer matrix is non-hermitian, 
    or force to do power projection 
"""
function non_hermitian_evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
                                init::Union{CMPS, Nothing} = nothing, 
                                max_pow_step::Integer = 100,
                                show_trace::Bool = false)
    """
        ResultFolder: classified by Model parameters(interaction, width), store CMPS and Obsv files
        CMPSResultFolder: CMPS information, classified by bond dimension
        ChkpFolder: Check points
    """
    CMPSResultFolder = @sprintf "%s/bondD_%02i_CMPS" ResultFolder bondD
    ChkpFolder = @sprintf "%s/CheckPoint_beta_%.2f" CMPSResultFolder β
    ChkpPsiFolder = @sprintf "%s/cmps_psi" ChkpFolder
    ChkpLpsiFolder = @sprintf "%s/cmps_Lpsi" ChkpFolder
    ChkpOptPsiFolder = @sprintf "%s/opt_psi" ChkpFolder
    isdir(ResultFolder) || mkdir(ResultFolder)
    isdir(CMPSResultFolder) || mkdir(CMPSResultFolder) 
    isdir(ChkpFolder) || mkdir(ChkpFolder)
    isdir(ChkpPsiFolder) || mkdir(ChkpPsiFolder)
    isdir(ChkpLpsiFolder) || mkdir(ChkpLpsiFolder)
    isdir(ChkpOptPsiFolder) || mkdir(ChkpOptPsiFolder)

    if m.Ut === nothing
        ChkpOptLpsiFolder = @sprintf "%s/opt_Lpsi" ChkpFolder
        isdir(ChkpOptLpsiFolder) || mkdir(ChkpOptLpsiFolder)
    end

    Tm = m.Tmatrix
    
    #initiate cmps
    pow_step = 0
    if init === nothing 
        ψr = init_cmps(bondD, Tm, β, show_trace = show_trace)
        m.Ut === nothing ? ψl = init_cmps(bondD, transpose(Tm), β, show_trace = show_trace) : ψl = m.Ut * ψr
        ψr = ψr |> diagQ
        ψl = ψl |> diagQ
    else 
        #incomplete code!
        ψ = init
        ψ = ψ |> diagQ
    end

    ChkpPsiFile = @sprintf "%s/step_%03i.hdf5" ChkpPsiFolder pow_step
    ChkpLpsiFile = @sprintf "%s/step_%03i.hdf5" ChkpLpsiFolder pow_step
    saveCMPS(ChkpPsiFile, ψr)
    saveCMPS(ChkpLpsiFile, ψl)

    ChkpEngFile = "$(ChkpFolder)/Obsv_FECvS.txt"
    open(ChkpEngFile,"w") do cfile
        F = free_energy(ψl, ψr, Tm, β)
        E = energy(ψl, ψr, Tm, β)
        Cv = specific_heat(ψl, ψr, Tm, β)
        S = β * (E - F)    
        write(cfile, "step      free_energy           energy              specific_heat            entropy      \n")
        write(cfile, "----  -------------------  --------------------   -------------------  -------------------\n")
        EngString = @sprintf "%3i   %.16f   %.16f   %.16f   %.16f \n" pow_step F E Cv S
        write(cfile, EngString)
    end

    PsiFidelityFile = "$(ChkpFolder)/Fidelity_psi.txt"
    open(PsiFidelityFile,"w") do cfile
        write(cfile, "step   fidelity_initial      fidelity_final  \n")
        write(cfile, "----  ------------------   ------------------\n")
    end

    if m.Ut === nothing
        LpsiFidelityFile = "$(ChkpFolder)/Fidelity_Lpsi.txt"
        open(LpsiFidelityFile,"w") do cfile
            write(cfile, "step   fidelity_initial      fidelity_final  \n")
            write(cfile, "----  ------------------   ------------------\n")
        end
    end

    while pow_step < max_pow_step
        pow_step += 1
        if show_trace println(@sprintf "Power Step = %03i" pow_step) end

        res = compress_cmps(Tm * ψr, bondD, β, show_trace = show_trace)
        ψr = res.ψ
        if m.Ut === nothing
            res2 = compress_cmps(transpose(Tm) * ψl, bondD, β, show_trace = show_trace)
            ψl = res2.ψ
        else
            ψl = m.Ut * ψr
        end

        open(ChkpEngFile,"a") do cfile
            F = free_energy(ψl, ψr, Tm, β)
            E = energy(ψl, ψr, Tm, β)
            Cv = specific_heat(ψl, ψr, Tm, β)
            S = β * (E - F)
            EngString = @sprintf "%3i   %.16f   %.16f   %.16f   %.16f \n" pow_step F E Cv S
            write(cfile, EngString)
        end

        ChkpPsiFile = @sprintf "%s/step_%03i.hdf5" ChkpPsiFolder pow_step
        ChkpLpsiFile = @sprintf "%s/step_%03i.hdf5" ChkpLpsiFolder pow_step
        saveCMPS(ChkpPsiFile, ψr)
        saveCMPS(ChkpLpsiFile, ψl)

        open(PsiFidelityFile,"a") do cfile
            FidelityString = @sprintf "%3i   %.16f   %.16f \n" pow_step res.fidelity_initial res.fidelity_final
            write(cfile, FidelityString)
        end

        OptPsitFile = @sprintf "%s/step_%03i.txt" ChkpOptPsiFolder pow_step
        open(OptPsitFile, "w") do file
            write(file, res.optim_result)
            if res.optim_result.trace !== nothing write(file, res.optim_result.trace) end
        end

        if m.Ut === nothing
            open(LpsiFidelityFile,"a") do cfile
                FidelityString = @sprintf "%3i   %.16f   %.16f \n" pow_step res2.fidelity_initial res2.fidelity_final
                write(cfile, FidelityString)
            end

            OptLpsiFile = @sprintf "%s/step_%03i.txt" ChkpOptLpsiFolder pow_step
            open(OptLpsiFile, "w") do file
                write(file, res2.optim_result)
                if res2.optim_result.trace !== nothing write(file, res2.optim_result.trace) end
            end
        end
    end
    return 
end

#end    # module evaluate