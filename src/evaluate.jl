#module evaluate
"""
    Initiate via boundary cMPS
"""
function init_cmps(bondD::Int64, model::PhysModel, β::Real)
    Tm = model.Tmatrix
    phy_dim = model.phy_dim

    ψ = CMPS(Tm.Q, Tm.R)
    pow_step, remainder = divrem(log(phy_dim, bondD), 1)

    for i = 1:Integer(pow_step - 1) ψ = Tm * ψ end
    if remainder != 0
        ψ0 = Tm * ψ
        ψ  = cmps_compress(ψ0, bondD, β)
    end
    return ψ
end


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
    isdir(ResultFolder) || mkdir(ResultFolder)
    isdir(CMPSResultFolder) || mkdir(CMPSResultFolder) 

    #initiate cmps
    init === nothing ? ψ = init_cmps(bondD, m, β) : ψ = init
    ψ = ψ |> diagQ

    loss() = free_energy(CMPS(diag(ψ.Q)|> diagm, ψ.R), m.Tmatrix, β)
    #FluxOptTools: Zygote.refresh() currently needed when defining new adjoints
    Zygote.refresh() 
    lossfun, gradfun, fg!, p0 = optfuns(loss, Params([ψ.Q, ψ.R]))
    opt = Optim.optimize(Optim.only_fg!(fg!), p0, LBFGS(),Optim.Options(iterations=10000, store_trace=true))
    
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
            ψ  = cmps_compress(ψ0, bondD, β, init = ψ)
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