@with_kw struct EvaluateOptions{Tinit<:Union{AbstractCMPS, Nothing}, 
                                Tprocessor<:Processor, 
                                Testimator<:EstimatorType,
                                Thermitian<:Union{Bool, Nothing},
                                Tcontinue<:Union{Bool, Integer},
                                Tcompress<:CompressOptions,
                                Toptim<:Optim.Options}
    init::Tinit = nothing
    processor::Tprocessor = CPU
    trace_estimator::Testimator 
    hermitian::Thermitian = nothing
    Continue::Tcontinue = false
    max_pow_step::Integer = 100
    group::Integer = 1
    tag::String = Dates.format(now(), "yyyy-mm-dd")
    show_trace::Bool = false
    store_trace::Bool = true
    compress_options::Tcompress = CompressOptions(trace_estimator = trace_estimator,
                        show_trace = show_trace,
                        store_trace = store_trace, 
                        processor = processor)
    #optim option for hermitian_evaluate
    optim_options::Toptim = Optim.Options(f_tol = eps(), 
                                  g_tol = 1.e-8,
                                  iterations = 10000, 
                                  store_trace = store_trace,
                                  show_trace = show_trace, 
                                  show_every = 10)
end


#"""
#    when `opt::EvaluateOptions`, update the processor in opt.trace_estimator according to `opt.processor`
#"""
#update_processor(opt::EvaluateOptions{Tinit, Tprocessor, Testimator,
#                    Thermitian, Tcontinue, Tcompress, Toptim}
#                    ) where {Tinit, Tprocessor, Testimator <: Nothing,
#                    Thermitian, Tcontinue, Tcompress, Toptim} = opt

#function update_processor(opt::EvaluateOptions{Tinit, Tprocessor, Testimator,
#                            Thermitian, Tcontinue, Tcompress, Toptim}
#                            ) where {Tinit, Tprocessor, Testimator <: TraceEstimator,
#                            Thermitian, Tcontinue, Tcompress, Toptim}
#    @unpack trace_estimator = opt
#    if trace_estimator.options.processor == opt.processor
#        return opt
#    else
#        new_options = FTLMOptions(trace_estimator.options, opt.processor)
#        new_estimator = TraceEstimator(trace_estimator.estimator, options=new_options)
#        new_compress_options = CompressOptions(opt.compress_options, 
#                                    processor = opt.processor, 
#                                    trace_estimator = new_estimator)
#        return EvaluateOptions(opt, 
#                                trace_estimator = new_estimator,
#                                compress_options = new_compress_options)
#
#    end
#        return 
#end