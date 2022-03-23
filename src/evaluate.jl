#module evaluate
"""
    Evaluate PhysModel m when the hermiticty of its transfer matrix is unknown or not specified
"""
function evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
                    init::Union{CMPS, Nothing} = nothing, 
                    max_pow_step::Integer = 50,
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

    #initiate cmps
    init === nothing ? 
        ψ = init_cmps(bondD, m, β, show_trace= show_trace) : ψ = init

    ψ = ψ |> diagQ
    loss() = free_energy(CMPS(diag(ψ.Q)|> diagm, ψ.R), m.Tmatrix, β)
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
    return ψ, dict
end


"""
    Evaluate PhysModel m when its transfer matrix is non-hermitian, 
    or force to do power projection 
"""
function non_hermitian_evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
                                init::Union{CMPS, Nothing} = nothing, 
                                max_pow_step::Integer = 50,
                                show_trace::Bool = false)
    """
        ResultFolder: classified by Model parameters(interaction, width), store CMPS and Obsv files
        CMPSResultFolder: CMPS information, classified by bond dimension
        ChkpFolder: Check points
    """
    CMPSResultFolder = @sprintf "%s/bondD_%02i_CMPS" ResultFolder bondD
    ChkpFolder = @sprintf "%s/CheckPoint_beta_%.2f" CMPSResultFolder β
    isdir(ResultFolder) || mkdir(ResultFolder)
    isdir(CMPSResultFolder) || mkdir(CMPSResultFolder) 
    isdir(ChkpFolder) || mkdir(ChkpFolder)

    Tm = m.Tmatrix

    #initiate cmps
    pow_step = 0
    init === nothing ? 
        ψ = init_cmps(bondD, m, β, show_trace = show_trace) : ψ = init
    ψ = ψ |> diagQ

    ψ0 = CMPS(Tm.Q, Tm.R)
    while size(ψ0.Q)[1] < bondD ψ0 = Tm * ψ0 end
    fidelity_initial = 1.0
    fidelity_final = fidelity(ψ, ψ0, β, Normalize = true)

    ChkpStateFile = @sprintf "%s/cmps_step_%03i.hdf5" ChkpFolder pow_step
    saveCMPS(ChkpStateFile, ψ)

    ChkpEngFile = "$(ChkpFolder)/Obsv_FECvS.txt"
    open(ChkpEngFile,"w") do cfile
        F = free_energy(m.Ut*ψ, ψ, Tm, β)
        E = energy(m.Ut*ψ, ψ, Tm, β)
        Cv = specific_heat(m.Ut*ψ, ψ, Tm, β)
        S = β * (E - F)    
        write(cfile, "step      free_energy          energy            specific_heat        entropy       \n")
        write(cfile, "---- ------------------- ------------------- ------------------- -------------------\n")
        EngString = @sprintf "%3i   %.16f   %.16f   %.16f   %.16f \n" pow_step F E Cv S
        write(cfile, EngString)
    end

    ChkpFidelityFile = "$(ChkpFolder)/Fidelity.txt"
    open(ChkpFidelityFile,"w") do cfile
        write(cfile, "step      free energy       fidelity_initial    fidelity_final\n")
        write(cfile, "------   ----------------   ----------------   ----------------\n")
        FidelityString = @sprintf "%3i   %.16f   %.16f \n" pow_step fidelity_initial fidelity_final
        write(cfile, FidelityString)
    end

    while pow_step < max_pow_step
        pow_step += 1
        res = compress_cmps(Tm * ψ, bondD, β, show_trace = show_trace)
        ψ = res.ψ

        ChkpStateFile = @sprintf "%s/cmps_step_%03i.hdf5" ChkpFolder pow_step
        saveCMPS(ChkpStateFile, ψ)

        open(ChkpEngFile,"a") do cfile
            F = free_energy(m.Ut*ψ, ψ, Tm, β)
            E = energy(m.Ut*ψ, ψ, Tm, β)
            Cv = specific_heat(m.Ut*ψ, ψ, Tm, β)
            S = β * (E - F)
            EngString = @sprintf "%3i   %.16f   %.16f   %.16f   %.16f \n" pow_step F E Cv S
            write(cfile, EngString)
        end

        open(ChkpFidelityFile,"w") do cfile
            FidelityString = @sprintf "%3i   %.16f   %.16f \n" pow_step res.fidelity_initial res.fidelity_final
            write(cfile, FidelityString)
        end

        OptResultFile = @sprintf "%s/cmps_opt_step_%03i.txt" ChkpFolder pow_step
        open(OptResultFile, "w") do file
            write(file, res.optim_result)
            if res.optim_result.trace !== nothing write(file, res.optim_result.trace) end
        end
    end

    # calculate thermal dynamic quanties
    dict = Dict()
    dict["F"] = free_energy(m.Ut*ψ, ψ, Tm, β)
    dict["E"] = energy(m.Ut*ψ, ψ, Tm, β)
    dict["Cv"] = specific_heat(m.Ut*ψ, ψ, Tm, β)
    dict["S"] = β * (dict["E"] - dict["F"])
    ResultFile = @sprintf "%s/beta_%.2f.hdf5" CMPSResultFolder β
    saveCMPS(ResultFile, ψ, dict)
    return ψ, dict
end

#end    # module evaluate