"""
    MeraUpdateOptions

Struct of options in adaptive MERA update function.
...
# Arguments
- `atol`: tolerance of absolute difference of the fidelity between two normalized CMPS and `1.0`, 
    defalt ``atol=1.e-5``.
- `loss_atol`: tolerance of absolute difference of the value of loss function between two MERA update steps,
    defalt ``loss_atol = 1.e-12``.
- `maxiter`: maximum iteration steps, defalt ``maxiter = 50``.
- `interpolate::Bool`: whether to interpolate between two MERA update steps,
    defalt ``interpolate = true``.
- `store_trace::Bool`: store_trace or not, defalt ``store_trace = false``.
- `show_trace::Bool`: show_trace or not, defalt ``show_trace = false``.
- `trace_estimator`: method to evaluate `logtrexp` functions, no defalt value.
...
"""
@with_kw struct MeraUpdateOptions{Ttol<:Real, Te}
    atol::Ttol = 1.e-5
    loss_atol::Ttol = 1.e-12
    maxiter::Int64 = 50
    interpolate::Bool = true
    store_trace::Bool = false
    show_trace::Bool = false
    trace_estimator::Te 
end


"""
    MeraUpdateStep

Data structure to store MERA update step informations.
...
# Arguments
- `id` : Step ID.
- `θ`  : interpolate angle
- `ldiff`: loss function difference between current step and the last step.
- `fidel`: fidelity between CMPS at current step and the origional CMPS.
...
"""
@with_kw struct MeraUpdateStep{Ti<:Integer, T<:Real, Tf<:Real}
    id::Ti
    θ::T
    ldiff::Tf
    fidel::Tf
end


"""
    MeraUpdateTrace = Vector{MeraUpdateStep}
"""
MeraUpdateTrace = Vector{MeraUpdateStep}


"""
    MeraUpdateResult

Struct for MERA update result.
...
# Arguments
- `cmps:` the result CMPS.
- `trace::MeraUpdateTrace`: trace of MERA update processure.
"""
@with_kw struct MeraUpdateResult{T<:AbstractCMPS}
    cmps::T
    trace::MeraUpdateTrace
end


"""
    println(step::MeraUpdateStep)
"""
function Base.println(step::MeraUpdateStep)
    str = @sprintf "%03i   %.16f   %.16e   %.16f \n" step.id step.θ step.ldiff step.fidel
    println(str)
    return
end


function _header(::MeraUpdateTrace)
    str = """
    -----------------------------MERA update------------------------------
    step           θ             loss difference           fidelity      
    ----  -------------------  --------------------   -------------------
    """
    return str
end

"""
    write(io::IO, mtrace::MeraUpdateTrace)
"""
function Base.write(io::IO, mtrace::MeraUpdateTrace)
    print(io, _header(mtrace))
    for state in mtrace
        writedlm(io, [state.id state.θ state.ldiff state.fidel])
    end
    return
end

