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
    tol::Float64 = eps(), maxiter::Integer=50, interpolate::Bool = false)
    _, v = symeigen(ψ0. Q)
    p_current = v[:, end-χ+1:end]
    ψ = project(ψ0, p_current)

    loss_previous = 9.9e9
    loss_diff = 9.9e9
    step = 0

    trace = MeraUpdateTrace()
    while step < maxiter
        θ = π
        
        loss() = logfidelity(project(ψ0, p_current), ψ0, β)
        loss_current = loss()
        loss_diff = abs(loss_current - loss_previous)
        loss_previous = loss_current

        #store_trace
        step_info = MeraUpdateStep(step, θ, loss_diff, fidelity(ψ, ψ0, β, Normalize = true))
        push!(trace, step_info)
        
        if loss_diff < tol break end

        #MERA update
        grad = Zygote.gradient(loss, Params([p_current]))[p_current]
        F = svd(grad)
        p_next = F.U * F.Vt
        step += 1
 
        #https://mathoverflow.net/questions/262560/natural-ways-of-interpolating-unitary-matrices
        #https://groups.google.com/forum/#!topic/manopttoolbox/2zhx67doXaU
        #interpolate between unitary matrices
        while interpolate
            θ = θ/2
            if θ < π/(1.9^12) #12-times bisection, cos(θ) = 0.9999989926433588
                interpolate = false
            else
                p_interpolate = interpolate_isometry(p_next, p_current, θ)
                loss_interpolate = logfidelity(project(ψ0, p_interpolate), ψ0, β)
                if loss_interpolate > loss_current
                    p_next = p_interpolate
                    interpolate = false
                end
            end     
        end
        p_current = p_next
        ψ = project(ψ0, p_current)
    end
    return ψ, trace
end


"""
    compress a CMPS(ψ0) to a target dimension(χ)
"""
function compress_cmps(ψ0::CMPS, χ::Integer, β::Real; 
    init::Union{CMPS, Nothing} = nothing, return_opt::Bool = false)
    χ0 = size(ψ0.Q)[1]
    length(size(ψ0.R)) == 2 ? vir_dim = 1 : vir_dim = size(ψ0.R)[3]
    if χ0 <= χ
        @warn "The bond dimension of the initial CMPS ≤ target bond dimension."
        res = (ψ, nothing)
    else
        init === nothing ? ψ = init_cmps(χ, vir_dim) : ψ = init
        if size(ψ.Q) != (χ, χ) 
            msg = "χ of init cmps should be $(χ) instead of $(size(ψ.Q)[1])"
            @error DimensionMismatch(msg)
        else
            ψ = ψ |> diagQ
            loss() = -logfidelity(CMPS(diag(ψ.Q)|> diagm, ψ.R), ψ0, β)
            p0, f, g! = optim_functions(loss, Params([ψ.Q, ψ.R]))
            opt = Optim.optimize(f, g!, p0, LBFGS(),Optim.Options(iterations=10000, store_trace = true))
            return_opt ? res = (ψ, opt) : res = (ψ, nothing)
        end
    end 
    return res
end


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
        ψ  = compress_cmps(ψ0, bondD, β)
    end
    return ψ
end
#end    # module cMPSCompress