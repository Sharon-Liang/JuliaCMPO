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
    return CMPS(Q,R)
end

"""
    Fidelity between the target cMPS |ψ⟩ and the origional cMPS T|r⟩:
        fidelity = ⟨ψ|T|r⟩/√(⟨ψ|ψ⟩)
        logfidelity(ψ, ψ0) = ln(Fd)
"""
logfidelity(ψ::CMPS, ψ0::CMPS, β::Real) = log_overlap(ψ, ψ0, β) - 0.5*log_overlap(ψ, ψ, β)

function fidelity(ψ::CMPS, ψ0::CMPS, β::Real; Normalize::Bool = false)
    if Normalize
        ψ = normalize(ψ, β)
        ψ0 = normalize(ψ0, β)
    end
    return  logfidelity(ψ, ψ0, β) |> exp
end


"""
    adaptive_mera_update: update the isometry using iterative SVD update with line search
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

function adaptive_mera_update(ψ0::CMPS, χ::Integer, β::Real; 
    options::MeraUpdateOptions = MeraUpdateOptions())
    step = 1
    logfidelity0 = logfidelity(ψ0, ψ0, β)
    loss(p_matrix) = logfidelity(project(ψ0, p_matrix), ψ0, β)

    _, v = symeigen(ψ0.Q)
    p_current = v[:, end-χ+1:end]
    
    loss_previous = 9.9e9
    loss_current = loss(p_current)
    adiff = abs(loss_current - logfidelity0)
    ldiff = abs(loss_current - loss_previous)
    
    trace = MeraUpdateTrace()
    step_info = MeraUpdateStep(step, π, ldiff, exp(-adiff))
    if options.store_trace push!(trace, step_info) end
    if options.show_trace println(step_info) end

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
function compress_cmps(ψ0::CMPS, χ::Integer, β::Real; 
    init::Union{CMPS, Nothing} = nothing, atol::Float64 = 1.e-5, 
    high_fidel_maxiter::Int64 = 1000, low_fidel_maxiter::Int64=5000,
    mera_update_options::MeraUpdateOptions = MeraUpdateOptions())
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
            abs(fidelity_initial - 1.0) < atol ? maxiter = high_fidel_maxiter : maxiter = low_fidel_maxiter

            ψ = ψ |> diagQ
            loss() = -logfidelity(CMPS(diag(ψ.Q)|> diagm, ψ.R), ψ0, β)
            p0, f, g! = optim_functions(loss, Params([ψ.Q, ψ.R]))

            optim_options = Optim.Options(iterations = maxiter,
                                show_trace = true, show_every = 10,
                                store_trace = true)
            optim_result = Optim.optimize(f, g!, p0, LBFGS(), optim_options)

            fidelity_final = fidelity(ψ, ψ0, β, Normalize = true)
        end
    end 
    return CompressResult(ψ, fidelity_initial, fidelity_final, optim_result)
end


"""
    Initiate via boundary cMPS
"""
function init_cmps(bondD::Int64, model::PhysModel, β::Real)
    Tm = model.Tmatrix
    ψ = CMPS(Tm.Q, Tm.R)
    while size(ψ.Q)[1] < bondD
        ψ = Tm * ψ
    end

    if size(ψ.Q)[1] > bondD 
        res = compress_cmps(ψ, bondD, β, low_fidel_maxiter=1000)
        ψ = res.ψ
    end
    return ψ
end
#end    # module cMPSCompress