#module evaluate
"""
    Evaluate PhysModel m when the hermiticty of its transfer matrix is unknown or not specified
"""
function evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
                    init::Union{CMPS, Nothing} = nothing, max_pow_step::Integer = 100,
                    hermitian::Union{Bool, Nothing} = nothing)
    hermitian === nothing ? hermitian = ishermitian(m.Tmatrix) : hermitian = hermitian
    if hermitian
        hermitian_evaluate(m, bondD, β, ResultFolder, init = init)
    else
        non_hermitian_evaluate(m, bondD, β, ResultFolder, init = init, max_pow_step = max_pow_step)
    end
end


"""
    Evaluate PhysModel m when its transfer matrix is hermitian
"""
function hermitian_evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; init::Union{CMPS, Nothing} = nothing)
    """
    ResultFolder: classified by Model parameters(interaction, width), store CMPS and Obsv files
    CMPSResultFolder: CMPS information, classified by bond dimension
    """
    CMPSResultFolder = @sprintf "%s/bondD_%02i_CMPS" ResultFolder bondD
    OptResultFolder = @sprintf "%s/bondD_%02i_Opt" CMPSResultFolder bondD
    isdir(ResultFolder) || mkdir(ResultFolder)
    isdir(CMPSResultFolder) || mkdir(CMPSResultFolder) 
    isdir(OptResultFolder) || mkdir(OptResultFolder)

    #initiate cmps
    init === nothing ? ψ = init_cmps(bondD, m, β) : ψ = init

    ψ = ψ |> diagQ
    loss() = free_energy(CMPS(diag(ψ.Q)|> diagm, ψ.R), m.Tmatrix, β)

    p0, f, g! = optim_functions(loss, Params([ψ.Q, ψ.R]))
    opt = Optim.optimize(f, g!, p0, LBFGS(),Optim.Options(iterations=10000))
    
    #save optimize result
    OptResultFile = @sprintf "%s/beta_%.2f.txt" OptResultFolder β
    open(OptResultFile, "w") do file
        write(file, opt)
    end

    # calculate thermal dynamic quanties
    dict = Dict()
    dict["F"] = minimum(opt)
    ResultFile = @sprintf "%s/beta_%.2f.hdf5" CMPSResultFolder β
    saveCMPS(ResultFile, ψ, dict)
    return ψ, dict
end


"""
    Evaluate PhysModel m when its transfer matrix is non-hermitian, 
    or force to do power projection 
"""
function non_hermitian_evaluate(m::PhysModel, bondD::Integer, β::Real, ResultFolder::String; 
                                init::Union{CMPS, Nothing} = nothing, max_pow_step::Integer = 100)
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

    #initiate cmps
    init === nothing ? ψ = init_cmps(bondD, m, β) : ψ = init
    ψ = ψ |> diagQ

    ChkpEngFile = "$(ChkpFolder)/Obsv_F.txt"
    open(ChkpEngFile,"w") do cfile
        for pow_step = 0 : max_pow_step
            ChkpStateFile = @sprintf "%s/cmps_step_%03i.hdf5" ChkpFolder pow_step
            saveCMPS(ChkpStateFile, ψ)
            writedlm(cfile, [pow_step free_energy(m.Ut*ψ, ψ, m.Tmatrix, β)])

            ψ0 = m.Tmatrix * ψ
            ψ, opt = compress_cmps(ψ0, bondD, β, init = ψ, return_opt = true)
            OptResultFile = @sprintf "%s/cmps_step_%03i.txt" ChkpFolder β
            open(OptResultFile, "w") do file
                write(file, opt)
            end
        end
    end

    # calculate thermal dynamic quanties
    dict = Dict()
    dict["F"] = free_energy(m.Ut*ψ, ψ, m.Tmatrix, β)

    ResultFile = @sprintf "%s/beta_%.2f.hdf5" CMPSResultFolder β
    saveCMPS(ResultFile, ψ, dict)
    return ψ, dict
end

#end    # module evaluate