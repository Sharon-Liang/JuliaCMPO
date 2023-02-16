"""
    CompressOptions

Structure of options when compressing a CMPS.
...
# Arguments
- `init`: initial guess of the compressed CMPS, one uses a randomly initiated CMPS when ``init=nothing``,
    defalt ``init = nothing``.
- `processor::Processor`: processor used in the calculation, defalt ``processor = CPU``.
- `store_trace::Bool`: store_trace or not, defalt ``store_trace = false``.
- `show_trace::Bool`: show_trace or not, defalt ``show_trace = false``.
- `trace_estimator`: method to evaluate `logtrexp` functions, no defalt value.
- `mera_update_options::MeraUpdateOptions`
- `optim_options::Optim.Options`
...
"""
@with_kw struct CompressOptions{Ti, Te}
    init::Ti = nothing
    processor::Processor = CPU
    show_trace::Bool = false
    store_trace::Bool = true
    trace_estimator::Te
    mera_update_options::MeraUpdateOptions = 
        MeraUpdateOptions(trace_estimator = trace_estimator,
                            show_trace = show_trace,
                            store_trace = store_trace)
    # The same as scipy L-BFGS-B
    optim_options::Optim.Options = 
        Optim.Options(f_tol = 2.220446049250313e-9, 
                      g_tol = 1.e-5,
                      iterations = 100,
                      store_trace = store_trace,
                      show_trace = show_trace, 
                      show_every = 10)
end


"""
    CompressResult

Structure to store the compressed results.
...
# Arguments
- `cmps`: the compressed CMPS.
- `ifidel`: the initial fidelity between the compressed CMPS and the origional CMPS before optimize.
- `ffidel`: the final fidelity between the compressed CMPS and the origional CMPS.
- `optim_result`: the result of `Optim.optim` function.
...
"""
@with_kw struct CompressResult{Ts<:AbstractCMPS, Tf<:Real, To}
    cmps::Ts
    ifidel::Tf
    ffidel::Tf
    optim_result::To
end