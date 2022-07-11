const EstimatorType = Union{TraceEstimator, Nothing}

"""
    `MeraUpdateOptions` struct:
    `atol`: tolerance of absolute difference of fidelity(ψ, ψ0, β, Normalize = true) and 1.0
    `ldiff_tol`: tolerance of absolute difference of the value of loss function between two MERA update steps
    `maxiter`: maximum iteration steps
    `interpolate`: interpolate between two steps of MERA update or nothing
    `store_trace`: store_trace or not 
    `show_trace`: show_trace or not
    `trace_estimator`: trace estimator
"""
@with_kw struct MeraUpdateOptions{Ttol<:Number, Tb<:Bool, 
                                  Testimator<:EstimatorType}
    atol::Ttol = 1.e-5
    ldiff_tol::Ttol = 1.e-12
    maxiter::Int64 = 100
    interpolate::Tb = true
    store_trace::Tb = false
    show_trace::Tb = false
    trace_estimator::Testimator = nothing
end


"""
    Data structure of MERA update steps
    `id` : Step ID
    `θ`  : interpolate angle
    `loss_diff`: loss function difference to last step
    `fidelity` : fidelity between `ψ` and the origional `ψ0`
"""
@with_kw struct MeraUpdateStep{Ti<:Integer, T<:Real, Tf<:Real}
    id::Ti
    θ::T
    loss_diff::Tf
    fidelity::Tf
end

"""
    `MeraUpdateTrace`
"""
MeraUpdateTrace = Vector{MeraUpdateStep}


"""
    `MeraUpdateResult`:
    `ψ`: updated CMPS
    `trace`: MeraUpdateTrace
"""
@with_kw struct MeraUpdateResult{T<:AbstractCMPS}
    ψ::T
    trace::MeraUpdateTrace
end


"""
    println MeraUpdateStep
"""
function Base.println(step::MeraUpdateStep)
    str = @sprintf "%03i   %.16f   %.16e   %.16f \n" step.id step.θ step.loss_diff step.fidelity
    println(str)
    return
end


"""
    write MERA update trace
"""
function Base.write(io::IO, tr::MeraUpdateTrace)
    @printf io "step          θ            loss_diff           fidelity    \n"
    @printf io "----  ----------------  ----------------   ----------------\n"
    for state in tr
        writedlm(io, [state.id state.θ state.loss_diff state.fidelity])
    end
    return
end

