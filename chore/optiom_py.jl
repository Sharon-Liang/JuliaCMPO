using cMPO
using Zygote, LinearAlgebra
using PyCall, SciPy; so = SciPy.optimize

function optim_functions_py(loss, pars::Zygote.Params)
    grads = Zygote.gradient(loss, pars)
    p0 = zeros(pars)
    copy!(p0, pars)

    gradient_function = function (w)
        copy!(pars, w)
        grads = Zygote.gradient(loss, pars)
        g = similar(w)
        copy!(g, grads)
    end

    loss_function = function (w)
        copy!(pars, w)
        loss()
    end
    return p0, loss_function, gradient_function
end

function compress_cmps_py(ψ0::CMPS, χ::Integer, β::Real; 
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
            println("loss function")

            p0, f, g! = optim_functions_py(loss, Params([ψ.Q, ψ.R]))
            println("gradient function")

            optim_options = Dict("disp"=> 10)
            optim_result = so.minimize(f, p0, method="L-BFGS-B", jac=g!, options = optim_options)

            fidelity_final = fidelity(ψ, ψ0, β, Normalize = true)
        end
    end 
    return CompressResult(ψ, fidelity_initial, fidelity_final, optim_result)
end


"""
    Initiate via boundary cMPS
"""

function init_cmps_py(bondD::Int64, model::PhysModel, β::Real)
    Tm = model.Tmatrix
    ψ = CMPS(Tm.Q, Tm.R)
    while size(ψ.Q)[1] < bondD
        ψ = Tm * ψ
    end

    if size(ψ.Q)[1] > bondD 
        res = compress_cmps_py(ψ, bondD, β, low_fidel_maxiter=1000)
        ψ = res.ψ
    end
    return ψ
end
