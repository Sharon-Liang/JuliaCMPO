@with_kw struct CompressOptions{Ti<:Union{AbstractCMPS, Nothing}, 
                         Tm<:MeraUpdateOptions,
                         To<:Optim.Options,
                         Ttrace<:EstimatorType}
    init::Ti = nothing
    show_trace::Bool = false
    store_trace::Bool = true
    trace_estimator::Ttrace = nothing
    mera_update_options::Tm = MeraUpdateOptions(MeraUpdateOptions(),
                            show_trace = show_trace,
                            store_trace = store_trace,
                            trace_estimator = trace_estimator)
    # The same as scipy L-BFGS-B
    optim_options::To = Optim.Options(f_tol = 2.220446049250313e-9, 
                                  g_tol = 1.e-5,
                                  iterations = 200,
                                  store_trace = store_trace,
                                  show_trace = show_trace, 
                                  show_every = 10)
end


"""
    compress a CMPS(ψ0) to a target dimension(χ)
"""
@with_kw struct CompressResult{T<:AbstractCMPS, Tf<:Real, To}
    ψ::T
    fidelity_initial::Tf
    fidelity_final::Tf
    optim_result::To
end