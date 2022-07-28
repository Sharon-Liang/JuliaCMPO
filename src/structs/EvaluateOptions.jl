"""
    EvaluateOptions

Struct of options in `evaluate` function.
...
# Arguments
- `init`: initial guess of the result CMPS, one uses a randomly initiated CMPS when ``init=nothing``,
    defalt ``init = nothing``. 
- `processor::Processor`: processor used in the calculation, defalt ``processor = CPU``.
- `store_trace::Bool`: store_trace or not, defalt ``store_trace = true``.
- `show_trace::Bool`: show_trace or not, defalt ``show_trace = false``.
- `trace_estimator`: method to evaluate `logtrexp` functions, no defalt value.
- `hermitian::Union{Bool, Nothing}`: use `hermitian_evaluate` to optimize free energy 
    or use `non_hermitian_evaluate` to apply power method. If `hermitian = nothing`, 
    the method is determined by the hermiticty of the CMPO. defalt `hermitian = nothing`.
- `tocontinue::Union{Bool, Nothing}`: Whether to continue calculation, or continue
    with the `tocontinue`-th step. defalt `tocontinue = false`.
- `max_pow_step`: maximum step of power iteration.
- `group`: Number of CMPOs to group together, normally it equals to the sublattices 
    in a unit cell, defalt `group = 1`
- `tag::String`: tag to the current calculation, defalt `taf = Dates.format(now(), "yyyy-mm-dd")`.
- `compress_options`: options to compress CMPSs.
- `optim_options`: options to optimize the free energy in `hermitian_evaluate` function.
...
"""
@with_kw struct EvaluateOptions{Ti,Te,Th<:Union{Bool, Nothing},Tc<:Union{Bool, Integer}}
    init::Ti = nothing
    processor::Processor = CPU
    show_trace::Bool = false
    store_trace::Bool = true
    trace_estimator::Te 
    hermitian::Th = nothing
    tocontinue::Tc = false
    max_pow_step::Integer = 100
    group::Integer = 1
    tag::String = Dates.format(now(), "yyyy-mm-dd")
    compress_options::CompressOptions = CompressOptions(trace_estimator = trace_estimator,
                        show_trace = show_trace,
                        store_trace = store_trace, 
                        processor = processor)
    #optim option for hermitian_evaluate
    optim_options::Optim.Options = Optim.Options(f_tol = eps(), 
                                  g_tol = 1.e-8,
                                  iterations = 2000, 
                                  store_trace = store_trace,
                                  show_trace = show_trace, 
                                  show_every = 10)
end
