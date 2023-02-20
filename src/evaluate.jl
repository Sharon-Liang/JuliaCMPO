#=
### *Hermitian* : *Variational Principle*
=#
"""
    variation_evaluate(Tₘ::CMPO, bondD::Int, β::Real[, init=nothing]; kwargs...)

Calculate the dominant eigen cMPS of a cMPO transfer matrix `Tₘ` at temperature `1/β` by variational principle, i.e. minimizing free energy.

Arguments:
* Tₘ    : the cMPO local tensor of the transfer matrix.
* bondD : target bond dimension of the result cMPS.
* β     : inverse temperature.
* init  : Union{Nothing, CMPS}, initial guess of the result cMPS.

Keyword Arguments:
* processor        : processor used.
* compress_options : CompressOptions used when initiating a cMPS
* optim_options    : Optim.Options used when variational minimizing free energy function.
"""
function variation_evaluate(Tₘ::CMPO, bondD::Integer, β::Real, init::Union{Nothing, CMPS} = nothing; 
    processor::Processor = CPU,
    compress_options::CompressOptions=CompressOptions(), 
    optim_options::Optim.Options = Optim.Options(f_tol = eps(), 
        g_tol = 1.e-8,
        iterations = 2000, 
        show_trace = true,
        show_every = 2)
    )

    solver = solver_function(processor)
    compress_options = CompressOptions(compress_options, processor=processor)

    Tₘ = solver(Tₘ)

    #Initiate cmps: initiate via boundary cMPS of init === noting
    init === nothing ? ψ = init_cmps(bondD, Tₘ, β; compress_options) : ψ = solver(init)
    ψ = diagQ(ψ)

    # loss function
    Qd = convert(Vector, diag(ψ.Q)) 
    R  = convert(Array, ψ.R)
    loss() = free_energy(solver(CMPS(diagm(Qd), R)), Tₘ , β)
    p₀, f, g! = optim_functions(loss, Zygote.Params([Qd, R]))

    println("Start optimizing: β = $(β)")
    res = Optim.optimize(f, g!, p₀, LBFGS(), optim_options)
    println(res)

    return CMPS(diagm(Qd), R) |> solver
end


"""
    variation_evaluate(Tₘ::CMPO, bondD::Int, βlist::Vector{<:Real}; kwargs...)

Calculate the dominant eigen cMPS of a cMPO transfer matrix `Tₘ` at a list of temperatures by variational principle, i.e. minimizing free energy. 

Note:
* Temperatures are manully sorted from high to low. 
* We use the result cMPS of the last step to initiating that of the next step.

Arguments:
* Tₘ    : the cMPO local tensor of the transfer matrix.
* bondD : target bond dimension of the result cMPS.
* βlist : inverse temperature list.
* init  : ``Union{Nothing, CMPS}``, initial guess of the result cMPS.

Keyword Arguments:
* processor        : processor used.
* compress_options : ``CompressOptions`` used when initiating a cMPS
* optim_options    : ``Optim.Options`` used when variational minimizing free energy function.
* obsv_functions   : ``Vector{<:Function}``, a list of thermaldynamic quanties are calculated.
* result_folder    : ``String``, folder to save the outputs, defalt = `"."` .
"""
function variation_evaluate(Tₘ::CMPO, bondD::Integer, βlist::Vector{<:Real}, init::Union{Nothing, CMPS} = nothing; 
    processor::Processor = CPU,
    compress_options::CompressOptions=CompressOptions(), 
    optim_options::Optim.Options = Optim.Options(f_tol = eps(), 
        g_tol = 1.e-8,
        iterations = 2000, 
        store_trace = true),
    obsv_functions::Vector{<:Function} = [free_energy, energy, specific_heat, entropy],
    result_folder::String = "."
    )

    # Evaluate from high T to low T
    sort!(βlist)
    solver = solver_function(processor)
    compress_options = CompressOptions(compress_options, processor=processor)
    isdir(result_folder) || mkdir(result_folder)

    Tₘ = solver(Tₘ)
    #save setups
    jldopen(result_folder*"/setups.jld", "w") do file
        addrequire(file, Optim)
        write(file, "Tₘ", Tₘ)
        write(file, "bondD", bondD)
        write(file, "init", init)
        write(file, "processor", processor)
        write(file, "compress_options", compress_options)
        write(file, "optim_options", optim_options)
    end

    #Initiate cmps: initiate via boundary cMPS of init === noting
    init === nothing ? ψ₀ = init_cmps(bondD, Tₘ, βlist[1]; compress_options) : ψ₀ = solver(init)

    #Array to store values of observables
    obsvs = zeros(lastindex(βlist), lastindex(obsv_functions)+1)
    obsvs[:, 1] = βlist
    for i in eachindex(βlist)
        β = βlist[i]
        ψ = variation_evaluate(Tₘ, bondD, β, ψ₀; processor, compress_options, optim_options)

        #Calculate thermal dynamic quanties
        for j in eachindex(obsv_functions)
            obsvs[i, j+1] = obsv_functions[j](ψ, Tₘ, β)
        end

        #save cMPS
        cmps_file = @sprintf "%s/cmps_beta_%f.jld" result_folder β
        jldopen(cmps_file, "w") do file
            write(file, "ψ",  ψ)
        end

        ψ₀ = ψ
    end

    #save observables
    open(result_folder*"/obsvs.txt", "w") do file
        write(file, @sprintf "%-20s" "β")
        for j in eachindex(obsv_functions)
            write(file, @sprintf "%-20s" obsv_functions[j])
        end
        write(file, "\n")
        writedlm(file, obsvs)
    end
end



#=
### *Non-Hermitian* : *Power method*
=#
"""
    power_evaluate(Tₘ::CMPO, bondD::Int, β::Real[, init=nothing]; kwargs...)

Calculate the dominant eigen cMPS of a cMPO transfer matrix `Tₘ` at temperature `1/β` by power method.

Arguments:
* Tₘ    : the cMPO local tensor of the transfer matrix.
* bondD : target bond dimension of the result cMPS.
* β     : inverse temperature.
* init  : ``Union{Nothing, Tuple}``, initial guess of the result cMPS, `init` shold be provided as `(ψl, ψr)`.

Keyword Arguments:
* processor        : processor used.
* compress_options : ``CompressOptions`` used when initiating and compressing a cMPS.
* max_pow_step     : maxinum power step, defalt `= 100`.
* result_folder    : ``String``, folder to save the outputs, defalt = `"."` .
* to_group         : ``Int``, if to group two lattice sites. `to_group ≤ 1` means not to group. 
* to_shift         : ``Float64``, if to shift the spectrum of `Tₘ`, `to_shift = 0` means not to shift.
"""
function  power_evaluate(Tₘ::CMPO, bondD::Integer, β::Real, init::Union{Nothing, Tuple}=nothing; 
    processor::Processor,
    compress_options::CompressOptions=CompressOptions(), 
    result_folder::String = ".",
    max_pow_step:: Int = 100, 
    to_group::Int = 0
    )

    solver = solver_function(processor)
    compress_options = CompressOptions(compress_options, processor=processor)
    isdir(result_folder) || mkdir(result_folder)
    
    Tₘ = solver(Tₘ)
    g = 1
    while g < to_group #enlarged matrix: phy_dim -> phy_dim^group
        Tₘ *= Tₘ
        g += 1
    end

    #initiate cmps
    pow_step = 0.
    if init === nothing
        ψr = init_cmps(bondD, Tₘ, β; compress_options) |> diagQ
        ψl = init_cmps(bondD, transpose(Tₘ), β; compress_options) |> diagQ
    else
        ψl, ψr = solver(init[1]), solver(init[2])
    end 

    #Save check points
    println("Start optimizing: β = $(β)")
    ckpt_file = @sprintf "%s/cmps_ckpt_beta_%f.jld" result_folder β
    jldopen(ckpt_file, "w") do file
        write(file, string(pow_step), (ψl, ψr))
    end

    while pow_step < max_pow_step
        pow_step += 1
        println(@sprintf "Power Step = %03i" pow_step)

        ψr = compress_cmps(Tₘ * ψr, bondD, β; compress_options)
        ψl = compress_cmps(transpose(Tₘ) * ψl, bondD, β; compress_options)

        #save checkpoints
        jldopen(ckpt_file, "r+") do file
            write(file, string(pow_step), (ψl, ψr))
        end
    end

    return ψl, ψr
end


"""
    power_evaluate(Tₘ::CMPO, bondD::Int, βlist::Vector{<:Real}[, init=nothing]; kwargs...)

Calculate the dominant eigen cMPS of a cMPO transfer matrix `Tₘ` at temperature `1/β` by power method.

Arguments:
* Tₘ    : the cMPO local tensor of the transfer matrix.
* bondD : target bond dimension of the result cMPS.
* βlist : inverse temperature list.
* init  : ``Vector``, initial guess of the result cMPS. The length of it should match the length of `βlist`, i.e. initiate values of each elements in `βlist` should be specified.

Keyword Arguments:
* processor        : processor used.
* compress_options : ``CompressOptions`` used when initiating and compressing a cMPS.
* max_pow_step     : maxinum power step, defalt `= 100`.
* result_folder    : ``String``, folder to save the outputs, defalt = `"."` .
* to_group         : ``Int``, if to group two lattice sites. `to_group ≤ 1` means not to group. 
* to_shift         : ``Float64``, if to shift the spectrum of `Tₘ`, `to_shift = 0` means not to shift.
* obsv_functions   : ``Vector{<:Function}``, a list of thermaldynamic quanties are calculated.
"""
function  power_evaluate(Tₘ::CMPO, bondD::Integer, βlist::Vector{<:Real}, init::AbstractVector=fill(nothing, lastindex(βlist)); 
    processor::Processor = GPU,
    compress_options::CompressOptions=CompressOptions(), 
    result_folder::String = ".",
    max_pow_step:: Int = 100, 
    to_group::Int = 0,
    obsv_functions::Vector{<:Function} = [free_energy, energy, specific_heat, entropy]
    )

    # Evaluate from high T to low T
    sort!(βlist)
    solver = solver_function(processor)
    compress_options = CompressOptions(compress_options, processor=processor)
    isdir(result_folder) || mkdir(result_folder)
    
    Tₘ = solver(Tₘ)
    g = 1
    while g < to_group #enlarged matrix: phy_dim -> phy_dim^group
        Tₘ *= Tₘ
        g += 1
    end

    #save setups
    jldopen(result_folder*"/setups.jld", "w") do file
        addrequire(file, Optim)
        write(file, "Tₘ", Tₘ)
        write(file, "bondD", bondD)
        write(file, "init", init)
        write(file, "processor", processor)
        write(file, "compress_options", compress_options)
        write(file, "max_pow_step", max_pow_step)
        write(file, "to_group", to_group)
    end

    #Array to store values of observables
    obsvs = zeros(lastindex(βlist), lastindex(obsv_functions)+1)
    obsvs[:,1] = βlist
    for i in eachindex(βlist)
        β = βlist[i]

        ψl, ψr = power_evaluate(Tₘ, bondD, β, init[i]; to_group = 0, processor, compress_options, result_folder, max_pow_step)
    
        #Calculate thermal dynamic quanties
        for j in eachindex(obsv_functions)
            obsvs[i, j+1] = obsv_functions[j](ψl, ψr, Tₘ, β)/g
        end
    
        #save cMPS
        cmps_file = @sprintf "%s/cmps_beta_%f.jld" result_folder β
        jldopen(cmps_file, "w") do file
            write(file, "ψlr",  (ψl, ψr))
        end
    end

    #save observables
    open(result_folder*"/obsvs.txt", "w") do file
        write(file, @sprintf "%-20s" "β")
        for j in eachindex(obsv_functions)
            write(file, @sprintf "%-20s" obsv_functions[j])
        end
        write(file, "\n")
        writedlm(file, obsvs)
    end
end
