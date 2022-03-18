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
init_cmps(χ::Int64) = init_cmps(χ, 1)

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
    atol::Float64 = 1.e-5, ldiff_tol::Float64 = 1.e-12, maxiter::Integer=50, 
    interpolate::Bool = true,
    store_trace::Bool = false, show_trace::Bool=false)
    """
        atol: tolerance of absolute difference of fidelity(ψ, ψ0, β, Normalize = true) and 1.0
        ldiff_tol: tolerance of absolute difference of the value of loss function between two MERA update steps
    """
    #Initiate: step = 0
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
    if store_trace
        step_info = MeraUpdateStep(step, π, 9.9e9, exp(-adiff))
        push!(trace, step_info)
    end
    if show_trace println(step_info) end

    while step < maxiter
        step += 1   
        grad = Zygote.gradient(x -> loss(x), p_current)[1]
        F = svd(grad)
        p_next = F.U * F.Vt
 
        #https://mathoverflow.net/questions/262560/natural-ways-of-interpolating-unitary-matrices
        #https://groups.google.com/forum/#!topic/manopttoolbox/2zhx67doXaU
        #interpolate between unitary matrices
        θ = π
        proceed = interpolate
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
        if store_trace
            step_info = MeraUpdateStep(step, θ, ldiff, exp(-adiff))
            push!(trace, step_info)
        end 
        if show_trace println(step_info) end

        if adiff < atol || ldiff < ldiff_tol break end
    end
    ψ = project(ψ0, p_current)
    store_trace ? (return ψ, trace) : return ψ
end


"""
    compress a CMPS(ψ0) to a target dimension(χ)
"""
function compress_cmps(ψ0::CMPS, χ::Integer, β::Real; 
    init::Union{CMPS, Nothing} = nothing, atol::Float64 = 1.e-5, 
    return_opt::Bool = false, show_trace::Bool = false)
    χ0 = size(ψ0.Q)[1]
    length(size(ψ0.R)) == 2 ? vir_dim = 1 : vir_dim = size(ψ0.R)[3]
    if χ0 <= χ
        @warn "The bond dimension of the initial CMPS ≤ target bond dimension."
        res = (ψ, nothing)
    else
        init === nothing ? ψ = adaptive_mera_update(ψ0, χ, β) : ψ = init
        if size(ψ.Q) != (χ, χ) 
            msg = "χ of init cmps should be $(χ) instead of $(size(ψ.Q)[1])"
            @error DimensionMismatch(msg)
        else
            if abs(fidelity(ψ, ψ0, β, Normalize = true) - 1.0) < atol
                res = (ψ, nothing)
            else
                ψ = ψ |> diagQ
                loss() = -logfidelity(CMPS(diag(ψ.Q)|> diagm, ψ.R), ψ0, β)
                p0, f, g! = optim_functions(loss, Params([ψ.Q, ψ.R]))
                opt = Optim.optimize(f, g!, p0, LBFGS(),
                    Optim.Options(iterations=3000, store_trace = true, 
                    show_trace = show_trace, show_every = 10))
                return_opt ? res = (ψ, opt) : res = (ψ, nothing)
            end
        end
    end 
    return res
end


"""
    Initiate via boundary cMPS
"""
function init_cmps(bondD::Int64, model::PhysModel, β::Real)
    Tm = model.Tmatrix
    ψ = CMPS(Tm.Q, Tm.R)
    while size(ψ.Q)[1] <= bondD
        ψ = Tm * ψ
    end

    if size(ψ.Q)[1] > bondD 
        ψ, _= compress_cmps(ψ, bondD, β)
    end
    return ψ
end
#end    # module cMPSCompress