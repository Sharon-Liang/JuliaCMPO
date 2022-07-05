#module cMPSCompress
"""
    Randomly initiate a cMPS with bond dimension χ and virtual bond dimension vD
"""
function init_cmps(χ::Int64, vD::Int64 = 1; 
    hermitian::Bool = true, dtype::DataType = Float64)
    Q = rand(dtype, χ, χ)
    if vD == 1
        R = rand(dtype, χ, χ)
    else
        R = rand(dtype, χ, χ, vD)
    end

    if hermitian
        Q = Q |> symmetrize
        for d = 1:vD
            R[:,:,d] = R[:,:,d] |> symmetrize
        end
    end
    return CMPS_generate(Q,R)
end


"""
    adaptive_mera_update: update the isometry using iterative SVD update with line search.
    `interpolate_isometry(p1,p2,θ)`: interpolate between two isometries
"""
function interpolate_isometry(p1::AbstractMatrix, p2::AbstractMatrix, θ::Real)
    """ interpolate two isometries
        θ = π/2 : mix = p1
        θ = 0   : mix = p2
    """
    mix = sin(θ) * p1 + cos(θ) * p2
    F = svd(mix)
    return F.U * F.Vt 
end

"""
    Data structure of MERA update steps
    `id` : Step ID
    `θ`  : interpolate angle
    `loss_diff`: loss function difference to last step
    `fidelity` : fidelity between `ψ` and the origional `ψ0`
"""
@with_kw struct MeraUpdateStep{Ti<:Integer, T<:Real, Tf<:Real}
    id::Ti
    θ::T
    loss_diff::Tf
    fidelity::Tf
end
MeraUpdateTrace = Vector{MeraUpdateStep}

@with_kw struct MeraUpdateResult{T<:AbstractCMPS}
    ψ::T
    trace::MeraUpdateTrace
end


"""
atol: tolerance of absolute difference of fidelity(ψ, ψ0, β, Normalize = true) and 1.0
ldiff_tol: tolerance of absolute difference of the value of loss function between two MERA update steps
"""
@with_kw struct MeraUpdateOptions{T1<:Number, T2<:Bool}
    atol::T1 = 1.e-5
    ldiff_tol::T1 = 1.e-12
    maxiter::Int64 = 100
    interpolate::T2 = true
    store_trace::T2 = false
    show_trace::T2 = false
end

function adaptive_mera_update(ψ0::AbstractCMPS, χ::Integer, β::Real; 
    options::MeraUpdateOptions = MeraUpdateOptions())
    step = 1
    #logfidelity0 = logfidelity(ψ0, ψ0, β)
    logfidelity0 = 9.9e9
    loss(p_matrix) = logfidelity(project(ψ0, p_matrix), ψ0, β)

    _, v = eigensolver(ψ0.Q |> symmetrize)
    p_current = v[:, end-χ+1:end]
    
    loss_previous = 9.9e9
    loss_current = loss(p_current)
    adiff = abs(loss_current - logfidelity0)
    ldiff = abs(loss_current - loss_previous)
    
    trace = MeraUpdateTrace()
    step_info = MeraUpdateStep(step, π, ldiff, exp(-adiff))
    if options.store_trace push!(trace, step_info) end
    if options.show_trace
        println("-----------------------------MERA update------------------------------")
        println("step           θ                 loss_diff             fidelity      \n")
        println("----  -------------------  --------------------   -------------------\n")
        println(step_info) 
    end

    while step < options.maxiter
        step += 1   
        grad = Zygote.gradient(x -> loss(x), p_current)[1]
        F = svd(grad)
        p_next = F.U * F.Vt
 
        #https://mathoverflow.net/questions/262560/natural-ways-of-interpolating-unitary-matrices
        #https://groups.google.com/forum/#!topic/manopttoolbox/2zhx67doXaU
        #interpolate between unitary matrices
        θ = π
        proceed = options.interpolate
        while proceed
            θ = θ/2
            if θ < π/(1.9^12) #12-times bisection, cos(θ) = 0.9999989926433588
                p_next = p_current
                proceed = false
            else
                p_interpolate = interpolate_isometry(p_next, p_current, θ)
                loss_interpolate = logfidelity(project(ψ0, p_interpolate), ψ0, β)
                if loss_interpolate > loss_current
                    p_next = p_interpolate
                    proceed = false
                end
            end     
        end
        p_current = p_next
        
        loss_current = loss(p_current)
        adiff = abs(loss_current - logfidelity0)
        ldiff = abs(loss_current - loss_previous)
        loss_previous = loss_current

        #store_trace
        step_info = MeraUpdateStep(step, θ, ldiff, exp(-adiff))
        if options.store_trace push!(trace, step_info) end 
        if options.show_trace println(step_info) end

        if adiff < options.atol || ldiff < options.ldiff_tol break end
    end
    ψ = project(ψ0, p_current)
    return MeraUpdateResult(ψ, trace)
end


"""
    compress a CMPS(ψ0) to a target dimension(χ)
"""
@with_kw struct CompressResult{T<:AbstractCMPS, Tf<:Real}
    ψ::T
    fidelity_initial::Tf
    fidelity_final::Tf
    optim_result
end

function compress_cmps(ψ0::T, χ::Integer, β::Real; 
    init::Union{AbstractCMPS, Nothing} = nothing,
    show_trace::Bool = false,
    mera_update_options::MeraUpdateOptions = MeraUpdateOptions()) where T<:AbstractCMPS
    mera_update_options = MeraUpdateOptions(mera_update_options, show_trace=show_trace)
    if show_trace 
        println("----------------------------Compress CMPS-----------------------------") 
    end

    χ0 = size(ψ0.Q)[1]
    length(size(ψ0.R)) == 2 ? vir_dim = 1 : vir_dim = size(ψ0.R)[3]
    optim_result = nothing
    fidelity_initial = 1.
    fidelity_final = 1.

    if χ0 <= χ
        ψ = ψ0
        @warn "The bond dimension of the initial CMPS ≤ target bond dimension."
    else
        init === nothing ? 
            ψ = adaptive_mera_update(ψ0, χ, β, options = mera_update_options).ψ : 
            ψ = init
        if size(ψ.Q) != (χ, χ) 
            msg = "χ of init cmps should be $(χ) instead of $(size(ψ.Q)[1])"
            @error DimensionMismatch(msg)
        else
            fidelity_initial = fidelity(ψ, ψ0, β, Normalize = true)

            ψ = ψ |> diagQ
            dQ = convert(Vector, diag(ψ.Q))
            R =  convert(Array, ψ.R)

            T <: CuCMPS ? solver = gpu_solver : solver = cpu_solver
            loss() = -logfidelity(solver(CMPS_generate, diagm(dQ), R), ψ0, β)
            p0, f, g! = optim_functions(loss, Params([dQ, R]))

            # The same as scipy L-BFGS-B
            optim_options = Optim.Options(f_tol = 2.220446049250313e-9, g_tol = 1.e-5,
                                iterations = 200,
                                store_trace = true,
                                show_trace = show_trace, show_every = 10)
            optim_result = Optim.optimize(f, g!, p0, LBFGS(), optim_options)
            
            ψ = solver(CMPS_generate, diagm(dQ), R)
            fidelity_final = fidelity(ψ, ψ0, β, Normalize = true)
        end
    end 
    return CompressResult(ψ, fidelity_initial, fidelity_final, optim_result)
end


"""
    Initiate via boundary cMPS
"""
function init_cmps(bondD::Integer, Tm::AbstractCMPO, β::Real; 
                    show_trace::Bool = false)
    if show_trace 
        println("----------------------------Initiate CMPS-----------------------------") 
    end
    ψ = CMPS_generate(Tm.Q, Tm.R)
    while size(ψ.Q)[1] < bondD
        ψ = Tm * ψ 
    end

    if size(ψ.Q)[1] > bondD 
        res = compress_cmps(ψ, bondD, β, show_trace=show_trace)
        ψ = res.ψ
    end
    return ψ
end

#end    # module cMPSCompress