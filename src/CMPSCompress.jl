"""
    `compress_cmps`:
"""
function compress_cmps(ψ0::AbstractCMPS{T, S, U}, χ::Integer, β::Real; 
    options::CompressOptions=CompressOptions(trace_estimator=nothing)) where {T,S,U}
    @unpack (init, show_trace, mera_update_options, 
             optim_options, trace_estimator, processor) = options
    if show_trace 
        println("----------------------------Compress CMPS-----------------------------") 
    end
    solver = solver_function(processor)
    
    χ0 = size(ψ0.Q, 1)
    U <: AbstractMatrix ? vir_dim = 1 : vir_dim = size(ψ0.R, 3)
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
            msg = "χ of init cmps should be $(χ) instead of $(size(ψ.Q, 1))"
            @error DimensionMismatch(msg)
        else
            fidelity_initial = fidelity(ψ, ψ0, β, trace_estimator, Normalize = true)

            ψ = ψ |> diagQ
            dQ = convert(Vector, diag(ψ.Q))
            R =  convert(Array, ψ.R)

            loss() = -logfidelity(solver(CMPS_generate, diagm(dQ), R), ψ0, β, trace_estimator)
            p0, f, g! = optim_functions(loss, Params([dQ, R]))

            optim_result = Optim.optimize(f, g!, p0, LBFGS(), optim_options)
            ψ = solver(CMPS_generate, diagm(dQ), R)
            fidelity_final = fidelity(ψ, ψ0, β, trace_estimator, Normalize = true)
        end
    end 
    return CompressResult(ψ, fidelity_initial, fidelity_final, optim_result)
end
