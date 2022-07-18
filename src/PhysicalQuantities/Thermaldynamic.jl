#module ThermaldynamicQuanties
"""
    make operator ⊢o⊣
"""
function make_operator(Op::AbstractMatrix{T}, dim::Int64) where T
    eye = Matrix{T}(I, dim, dim)
    return eye ⊗ Op ⊗ eye
end

function make_operator(Op::AbstractMatrix{T}, ψ::AbstractCMPS) where T
    eye = Matrix{T}(I, size(ψ.Q))
    return eye ⊗ Op ⊗ eye
end


"""
    The thermal average of local opeartors ⊢o⊣ with respect to K = ψl * W * ψr
"""
function thermal_average(Op::AbstractMatrix, ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real)
    K = ψl * W * ψr 
    e, v = eigensolver(K)
    e = -β*e
    m = maximum(e); e = e .- m
    Op = v' * Op * v
    den = exp.(e) |> sum
    num = exp.(e) .* diag(Op) |> sum
    return num/den   
end
thermal_average(Op::AbstractMatrix, ψ::AbstractCMPS, W::AbstractCMPO, β::Real) = thermal_average(Op, ψ, ψ, W, β)



"""
    The thermal average of local opeartors ⊢o⊣ with respect to K = ψ * ψ
"""
function thermal_average(Op::AbstractMatrix, ψl::AbstractCMPS, ψr::AbstractCMPS, β::Real)
    K = ψl * ψr 
    e, v = eigensolver(K)
    e = -β*e
    m = maximum(e) ; e = e .- m
    Op = v' * Op * v
    den = exp.(e) |> sum
    num = exp.(e) .* diag(Op) |> sum
    return num/den
end
thermal_average(Op::AbstractMatrix, ψ::AbstractCMPS, β::Real) =thermal_average(Op, ψ, ψ, β)


"""
    Free energy: F = -1/(βL) lnZ
"""
function free_energy(ψl::AbstractCMPS, ψr::AbstractCMPS, 
                     W::AbstractCMPO, β::Real, trace_estimator::EstimatorType = nothing)
    res = log_overlap(ψl, W * ψr, β, trace_estimator) - log_overlap(ψl, ψr, β, trace_estimator)
    return -res/β
end
free_energy(ψ::AbstractCMPS, W::AbstractCMPO, β::Real, trace_estimator::EstimatorType = nothing) = 
    free_energy(ψ, ψ, W, β, trace_estimator)


"""
    Energy density: E = -∂lnZ/∂β
"""
function energy(ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real)
    K = ψl * W * ψr 
    H = ψl * ψr 
    res = thermal_average(K, ψl, ψr, W, β) - thermal_average(H, ψl, ψr, β)
    return res
end
energy(ψ::AbstractCMPS, W::AbstractCMPO, β::Real) = energy(ψ, ψ, W, β)


"""
    Specific heat: Cv = -β^2 ∂E/∂β
"""
function specific_heat(ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real; 
                        method::Symbol = :adiff)
    if method == :adiff
        K = ψl * W * ψr 
        H = ψl * ψr 
        K2 = K * K
        H2 = H * H
        c = thermal_average(K2, ψl, ψr, W, β) - thermal_average(K, ψl, ψr, W, β)^2
        c -= thermal_average(H2, ψl, ψr, β) - thermal_average(H, ψl, ψr, β)^2
    elseif method == :ndiff
        minus_e = b -> -energy(ψl, ψr, W, b)
        c = central_fdm(5, 1)(minus_e, β)
    else @error "method should be :adiff or :ndiff"
    end
    return β^2 * c
end
specific_heat(ψ::AbstractCMPS, W::AbstractCMPO, β::Real; 
                method::Symbol = :adiff) = specific_heat(ψ, ψ, W, β, method = method)


"""
    Entropy: S = β×(E - F)
"""
function entropy(ψl::AbstractCMPS, ψr::AbstractCMPS, W::AbstractCMPO, β::Real)
    s = energy(ψl, ψr, W, β) - free_energy(ψl, ψr, W, β)
    return β*s
end
entropy(ψ::AbstractCMPS, W::AbstractCMPO, β::Real) = entropy(ψ, ψ, W, β)

#end  # module ThermaldynamicQuanties
