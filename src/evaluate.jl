#module evaluate

"""
    Evaluate PhysModel `m` when the hermiticty of its transfer matrix is unknown or not specified
"""
function evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
                    init = nothing, 
                    max_pow_step::Integer = 100,
                    hermitian::Union{Bool, Nothing} = nothing,
                    group::Integer = 1,
                    show_trace::Bool = false,
                    Continue::Union{Bool, Integer} = false,
                    tag = Dates.format(now(), "yyyy-mm-dd"),
                    solver::Function = cpu_solver)
    hermitian === nothing ? hermitian = ishermitian(m.Tmatrix) : hermitian = hermitian
    if hermitian
        hermitian_evaluate(m, bondD, β, ResultFolder, 
            init = init, 
            show_trace = show_trace,
            tag = tag, 
            solver = solver)
    else
        non_hermitian_evaluate(m, bondD, β, ResultFolder, 
            init = init, 
            max_pow_step = max_pow_step,
            show_trace = show_trace,
            Continue = Continue,
            group = group,
            tag = tag, 
            solver = solver)
    end
end


"""
    Evaluate PhysModel m when its transfer matrix is hermitian
"""
function hermitian_evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
            init::Union{CMPS, Nothing} = nothing,
            show_trace::Bool=false,
            tag = Dates.format(now(), "yyyy-mm-dd"),
            solver::Function = cpu_solver)
    """
    ResultFolder: classified by Model parameters(interaction, width), store CMPS and Obsv files
    CMPSResultFolder: CMPS information, classified by bond dimension
    """
    CMPSResultFolder = @sprintf "%s/bondD_%02i_CMPS_%s" ResultFolder bondD tag
    OptResultFolder = @sprintf "%s/bondD_%02i_Opt_%s" ResultFolder bondD tag
    isdir(ResultFolder) || mkdir(ResultFolder)
    isdir(CMPSResultFolder) || mkdir(CMPSResultFolder) 
    isdir(OptResultFolder) || mkdir(OptResultFolder)

    Tm = solver(x->x, m.Tmatrix)
    #initiate cmps
    init === nothing ? 
        ψ = solver(x->init_cmps(bondD, x, β, show_trace = show_trace), Tm) : ψ = solver(x->x, init)

    ψ = solver(diagQ, ψ)

    dQ = convert(Vector, diag(ψ.Q))
    R  = convert(Array, ψ.R)

    loss() = free_energy(solver(CMPS_generate, consist_diagm(dQ), R), Tm, β)
    F_initial = loss()

    p0, f, g! = optim_functions(loss, Params([dQ, R]))
    optim_options = Optim.Options(f_tol = eps(), g_tol = 1.e-8,
                            iterations = 10000,
                            store_trace = true,
                            show_trace = show_trace, show_every = 10)
    optim_result = Optim.optimize(f, g!, p0, LBFGS(), optim_options)
    F_final = minimum(optim_result)

    ψ = solver(CMPS_generate, consist_diagm(dQ), R)
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
    saveCMPS(ResultFile, CTensor(ψ), dict)
    return  ψ, dict
end


"""
    Evaluate PhysModel m when its transfer matrix is non-hermitian, 
    or force to do power projection 
"""
function non_hermitian_evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
                                init = nothing, 
                                group::Integer = 1,
                                max_pow_step::Integer = 100,
                                show_trace::Bool = false,
                                Continue::Union{Bool, Integer} = false,
                                tag = Dates.format(now(), "yyyy-mm-dd"),
                                solver::Function = cpu_solver)
    """
        ResultFolder: classified by Model parameters(interaction, width), store CMPS and Obsv files
        CMPSResultFolder: CMPS information, classified by bond dimension
        ChkpFolder: Check points
    """
    Tm = solver(x->x, m.Tmatrix)
    g = 1
    while g < group #enlarged Tmatrix: phy_dim -> phy_dim^group
        Tm = solver(x-> x * Tm, m.Tmatrix)
        g+=1
    end

    CMPSResultFolder = @sprintf "%s/bondD_%02i_CMPS_%s" ResultFolder bondD tag
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

    if Continue == false || Continue == 0
        #CLEAR ChkpFolder
        if length(readdir(ChkpFolder)) != 0
            for file in readdir(ChkpFolder) rm("$(ChkpFolder)/$(file)",recursive=true) end
        end
        isdir(ChkpPsiFolder) || mkdir(ChkpPsiFolder)
        isdir(ChkpLpsiFolder) || mkdir(ChkpLpsiFolder)
        isdir(ChkpOptPsiFolder) || mkdir(ChkpOptPsiFolder)

        if m.Ut === nothing
            isdir(ChkpOptLpsiFolder) || mkdir(ChkpOptLpsiFolder)
        end

        #initiate cmps
        pow_step = 0
        if init === nothing 
            ψr = solver(x->init_cmps(bondD, x, β, show_trace = show_trace), Tm)
            m.Ut === nothing ? 
                ψl = solver(x->init_cmps(bondD, x, β, show_trace = show_trace), transpose(Tm)) : 
                ψl = solver(x -> x * ψr, m.Ut)
        else 
            ψr, ψl = solver(x->x, init)
        end

        ψr = solver(diagQ, ψr)
        ψl = solver(diagQ, ψl)
        ChkpPsiFile = @sprintf "%s/step_%03i.hdf5" ChkpPsiFolder pow_step
        ChkpLpsiFile = @sprintf "%s/step_%03i.hdf5" ChkpLpsiFolder pow_step
        saveCMPS(ChkpPsiFile, CTensor(ψr))
        saveCMPS(ChkpLpsiFile, CTensor(ψl))

        open(ChkpEngFile,"w") do cfile
            F = free_energy(ψl, ψr, Tm, β)/group
            E = energy(ψl, ψr, Tm, β)/group
            Cv = specific_heat(ψl, ψr, Tm, β)/group
            S = β * (E - F)    
            write(cfile, "step    free_energy/site        energy/site       specific_heat/site      entropy/site    \n")
            write(cfile, "----  -------------------  --------------------   -------------------  -------------------\n")
            EngString = @sprintf "%3i   %.16f   %.16f   %.16f   %.16f \n" pow_step F E Cv S
            write(cfile, EngString)
        end
    
        open(PsiFidelityFile,"w") do cfile
            write(cfile, "step   fidelity_initial      fidelity_final  \n")
            write(cfile, "----  ------------------   ------------------\n")
        end
    
        if m.Ut === nothing
            open(LpsiFidelityFile,"w") do cfile
                write(cfile, "step   fidelity_initial      fidelity_final  \n")
                write(cfile, "----  ------------------   ------------------\n")
            end
        end
    else
        if Continue == true
            pow_step = readdlm(ChkpEngFile, skipstart=2)[end,1]
        else
            pow_step = Continue
            EngData = readlines(ChkpEngFile, keep=true)
            open(ChkpEngFile, "w") do file
                for i = 1:pow_step+3 write(file, EngData[i]) end
            end
            PsiFidelity = readlines(PsiFidelityFile, keep=true)
            open(PsiFidelityFile, "w") do file
                for i = 1:pow_step+2 write(file, PsiFidelity[i]) end
            end
            if m.Ut === nothing
                LpsiFidelity = readlines(LpsiFidelityFile, keep=true)
                open(LpsiFidelityFile, "w") do file
                    for i = 1:pow_step+2 write(file, LpsiFidelity[i]) end
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

        res = solver(x ->compress_cmps(Tm * x, bondD, β, show_trace = show_trace), ψr)
        ψr = solver(x->x, res.ψ)
        if m.Ut === nothing
            res2 = solver(x->compress_cmps(transpose(Tm) * x, bondD, β, show_trace = show_trace), ψl)
            ψl = solver(x->x, res2.ψ)
        else
            ψl = solver(x -> x * ψr, m.Ut)
        end

        open(ChkpEngFile,"a") do cfile
            F = free_energy(ψl, ψr, Tm, β)/group
            E = energy(ψl, ψr, Tm, β)/group
            Cv = specific_heat(ψl, ψr, Tm, β)/group
            S = β * (E - F)
            EngString = @sprintf "%3i   %.16f   %.16f   %.16f   %.16f \n" pow_step F E Cv S
            write(cfile, EngString)
        end

        ChkpPsiFile = @sprintf "%s/step_%03i.hdf5" ChkpPsiFolder pow_step
        ChkpLpsiFile = @sprintf "%s/step_%03i.hdf5" ChkpLpsiFolder pow_step
        saveCMPS(ChkpPsiFile, CTensor(ψr))
        saveCMPS(ChkpLpsiFile, CTensor(ψl))

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