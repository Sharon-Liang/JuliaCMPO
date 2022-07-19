@with_kw struct CompressOptions{Ti<:Union{AbstractCMPS, Nothing}, 
                         Tm<:MeraUpdateOptions,
                         To<:Optim.Options,
                         Testimator<:EstimatorType,
                         Tprocessor<:Processor}
    init::Ti = nothing
    show_trace::Bool = false
    store_trace::Bool = true
    trace_estimator::Testimator
    processor::Tprocessor = CPU
    mera_update_options::Tm = MeraUpdateOptions(trace_estimator = trace_estimator,
                            show_trace = show_trace,
                            store_trace = store_trace
                            )
    # The same as scipy L-BFGS-B
    optim_options::To = Optim.Options(f_tol = 2.220446049250313e-9, 
                                  g_tol = 1.e-5,
                                  iterations = 200,
                                  store_trace = store_trace,
                                  show_trace = show_trace, 
                                  show_every = 10)
end

#"""
#    when `opt::CompressOptions`, update the processor in opt.trace_estimator according to `opt.processor`
#"""
#update_processor(opt::CompressOptions{Ti, Tm, To, Testimator, Tprocessor}
#    ) where {Ti, Tm, To, Testimator <: Nothing, Tprocessor} = opt

#function update_processor(opt::CompressOptions{Ti, Tm, To, Testimator, Tprocessor}
#                         ) where {Ti, Tm, To, Testimator <: TraceEstimator, Tprocessor}
#    @unpack trace_estimator = opt
#    if trace_estimator.options.processor == opt.processor
#        return opt
#    else
#        new_options = FTLMOptions(trace_estimator.options, opt.processor)
#        new_estimator = TraceEstimator(trace_estimator.estimator, options=new_options)
#        return CompressOptions(opt, trace_estimator = new_estimator)
#    end
#end

"""
    compress a CMPS(ψ0) to a target dimension(χ)
"""
@with_kw struct CompressResult{T<:AbstractCMPS, Tf<:Real, To}
    ψ::T
    fidelity_initial::Tf
    fidelity_final::Tf
    optim_result::To
end