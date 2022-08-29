"""
    compress_cmps(ψ0, χ, β; compress_options::CompressOptions)

Compress a CMPS `ψ0` to bond dimension `χ` at temperature `β`.
"""
function compress_cmps(ψ0::AbstractCMPS, χ::Integer, β::Real; compress_options::CompressOptions)
    @unpack (init, processor, show_trace, trace_estimator,
             mera_update_options, optim_options) = compress_options
    if show_trace println(_header("Compress CMPS")) end
    solver = solver_function(processor)
    
    optim_result = nothing
    ifidel, ffidel, ifnormalize = 9.9e9, 9.9e9, true

    if size(ψ0.Q, 1) <= χ
        ψ = ψ0
        @warn "The bond dimension of the initial CMPS ≤ target bond dimension."
    else
        init === nothing ? 
            ψ = adaptive_mera_update(ψ0, χ, β; mera_update_options).cmps : 
            ψ = init
        if size(ψ.Q) != (χ, χ) 
            msg = "χ of initial CMPS should be $(χ) instead of $(size(ψ.Q, 1))"
            @error DimensionMismatch(msg)
        else
            ifidel = fidelity(ψ, ψ0, β, trace_estimator; ifnormalize)

            ψ = ψ |> diagQ
            dQ = convert(Vector, diag(ψ.Q))
            R =  convert(Array, ψ.R)

            loss() = -logfidelity(solver(cmps_generate(diagm(dQ), R)), ψ0, β, trace_estimator)
            p0, f, g! = optim_functions(loss, Params([dQ, R]))

            optim_result = Optim.optimize(f, g!, p0, LBFGS(), optim_options)
            ψ = solver(cmps_generate, diagm(dQ), R)
            ffidel = fidelity(ψ, ψ0, β, trace_estimator; ifnormalize)
        end
    end 
    return CompressResult(ψ, ifidel, ffidel, optim_result)
end
