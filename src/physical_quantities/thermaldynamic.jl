#module ThermaldynamicQuanties

#TODO: wrong thermal_average, energy and entropy function
"""
    make_operator(Op, dim::Integer)
    make_operator(Op, ψ::Integer)

Generate a operator ``⊢o⊣``.
"""
function make_operator(Op, dim::Integer)
    eye = Matrix{eltype(Op)}(I, dim, dim)
    return eye ⊗ Op ⊗ eye
end

function make_operator(Op, ψ::AbstractCMPS)
    eye = oneunit(ψ.Q)
    return eye ⊗ Op ⊗ eye
end


""" 
    free_energy(ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β[, ::Nothing])
    free_energy(ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β, trace_estimator::TraceEstimator)

Calculate the total free_energy ``F = -\\frac{1}{β}lnZ``,
where ``β = 1/T`` is the inverse temperature.
"""
function free_energy(ψl::AbstractCMPS, ψr::AbstractCMPS, 
                     W::AbstractCMPO, β::Real, trace_estimator = nothing)
    res = log_overlap(ψl, W * ψr, β, trace_estimator) - log_overlap(ψl, ψr, β, trace_estimator)
    return -res/β
end
free_energy(ψ::AbstractCMPS, W::AbstractCMPO, β::Real, trace_estimator = nothing) = 
    free_energy(ψ, ψ, W, β, trace_estimator)


"""
    thermal_average(Op, ψl, ψr, W, β[, ::Nothing])
    thermal_average(Op, ψl, ψr, W, β, trace_estimator::TraceEstimator)

The thermal average of local opeartors ``⊢o⊣`` with respect to ``K = ψl * W * ψr``,
```math
    ⟨Op⟩ = \frac{Tr(Op e^{-βK})}{Z}
```
"""
thermal_average(Op, ψl, ψr, W, β) = thermal_average(Op, ψl, ψr, W, β, nothing)
thermal_average(Op, ψ::AbstractCMPS, W::AbstractCMPO, β::Real, trace_estimator = nothing) = thermal_average(Op, ψ, ψ, W, β, trace_estimator)

function thermal_average(Op::AbstractArray, ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real, ::Nothing)
    K = ψl * W * ψr |> tomatrix |> symmetrize
    e, v = eigensolver(K)
    e0 = minimum(e)
    e = e .- e0
    Λ = exp.(-β .* e)
    ave = ein"ik, kl, li -> i"(conj(v), Op, v)

    den = sum(Λ)
    num = reduce(+, map(*, Λ, ave))
    return num/den   
end
thermal_average(Op::CMPSMatrix, ψl, ψr, W, β, trace_estimator::Nothing) = thermal_average(tomatrix(Op), ψl, ψr, W, β, trace_estimator)

function thermal_average(Op, ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real, trace_estimator::TraceEstimator)
    K = CMPSMatrix(ψl, W * ψr)
    return FiniteTLanczos.thermal_average(K, β, Op; trace_estimator)
end


"""
    thermal_average(Op, ψl, ψr, β[, ::Nothing])
    thermal_average(Op, ψl, ψr, β, trace_estimator::TraceEstimator)

    The thermal average of local opeartors ``⊢o⊣`` with respect to ``K = ψ * ψ``
"""
thermal_average(Op, ψl, ψr, β) = thermal_average(Op, ψl, ψr, β, nothing)
thermal_average(Op, ψ::AbstractCMPS, β::Real, trace_estimator = nothing) = thermal_average(Op, ψ, ψ, β, trace_estimator)

#TODO: wrong energy function: ave wrong!
function thermal_average(Op::AbstractArray, ψl::AbstractCMPS, ψr::AbstractCMPS, β::Real, ::Nothing)
    K = ψl * ψr |> tomatrix |> symmetrize
    e, v = eigensolver(K)
    e0 = minimum(e)
    e = e .- e0
    Λ = exp.(-β .* e)
    ave = ein"ik, kl, li -> i"(conj(v), Op, v)

    den = sum(Λ)
    num = reduce(+, map(*, Λ, ave))
    return num/den  
end
thermal_average(Op::CMPSMatrix, ψl, ψr, β, trace_estimator::Nothing) = thermal_average(tomatrix(Op), ψl, ψr, β, trace_estimator)

function thermal_average(Op, ψl::AbstractCMPS, ψr::AbstractCMPS, β::Real, trace_estimator::TraceEstimator)
    K = CMPSMatrix(ψl, ψr)
    return FiniteTLanczos.thermal_average(K, β, Op; trace_estimator)
end


"""
    Energy density: E = -∂lnZ/∂β
"""
function energy(ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real, trace_estimator = nothing)
    K = ψl * (W * ψr) 
    H = ψl * ψr 
    res = thermal_average(K, ψl, ψr, W, β, trace_estimator) - thermal_average(H, ψl, ψr, β, trace_estimator)
    return res
end
energy(ψ::AbstractCMPS, W::AbstractCMPO, β::Real, trace_estimator=nothing) = energy(ψ, ψ, W, β, trace_estimator)


#"""
#    Specific heat: Cv = -β^2 ∂E/∂β
#"""
#function specific_heat(ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real; 
#                        method::Symbol = :adiff)
#    if method == :adiff
#        K = ψl * (W * ψr) 
#        H = ψl * ψr 
#        K2 = K * K
#        H2 = H * H
#        c = thermal_average(K2, ψl, ψr, W, β) - thermal_average(K, ψl, ψr, W, β)^2
#        c -= thermal_average(H2, ψl, ψr, β) - thermal_average(H, ψl, ψr, β)^2
#    elseif method == :ndiff
#        minus_e = b -> -energy(ψl, ψr, W, b)
#        c = central_fdm(5, 1)(minus_e, β)
#    else @error "method should be :adiff or :ndiff"
#    end
#    return β^2 * c
#end
#specific_heat(ψ::AbstractCMPS, W::AbstractCMPO, β::Real; 
#                method::Symbol = :adiff) = specific_heat(ψ, ψ, W, β, method = method)


"""
    Entropy: S = β×(E - F)
"""
function entropy(ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real, trace_estimator = nothing)
    s = energy(ψl, ψr, W, β, trace_estimator) - free_energy(ψl, ψr, W, β, trace_estimator)
    return β*s
end
entropy(ψ::AbstractCMPS, W::AbstractCMPO, β::Real, trace_estimator = nothing) = entropy(ψ, ψ, W, β, trace_estimator)

#end  # module ThermaldynamicQuanties
