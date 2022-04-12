#module ThermaldynamicQuanties
"""
    make operator ⊢o⊣
"""
function make_operator(Op::AbstractMatrix{T}, dim::Int64) where T
    eye = Matrix{T}(I, dim, dim)
    return eye ⊗ Op ⊗ eye
end

function make_operator(Op::AbstractMatrix{T}, ψ::CMPS) where T
    eye = Matrix{T}(I, size(ψ.Q))
    return eye ⊗ Op ⊗ eye
end


"""
    The thermal average of local opeartors ⊢o⊣ with respect to K = ψl * W * ψr
"""
function thermal_average(Op::AbstractMatrix, ψl::CMPS, ψr::CMPS, W::CMPO, β::Real;
                            device=:cpu)
    K = ψl * W * ψr 
    e, v = symeigen(-β*K, device=device)
    m = maximum(e); e = e .- m
    Op = v' * Op * v
    den = exp.(e) |> sum
    num = exp.(e) .* diag(Op) |> sum
    return num/den   
end
thermal_average(Op::AbstractMatrix, ψ::CMPS, W::CMPO, β::Real; device=:cpu) = thermal_average(Op, ψ, ψ, W, β, device=device)



"""
    The thermal average of local opeartors ⊢o⊣ with respect to K = ψ * ψ
"""
function thermal_average(Op::AbstractMatrix, ψl::CMPS, ψr::CMPS, β::Real; device =:cpu)
    K = ψl * ψr 
    e, v = symeigen(-β*K, device=device)
    m = maximum(e) ; e = e .- m
    Op = v' * Op * v
    den = exp.(e) |> sum
    num = exp.(e) .* diag(Op) |> sum
    return num/den
end
thermal_average(Op::AbstractMatrix, ψ::CMPS, β::Real; device=:cpu) =thermal_average(Op, ψ, ψ, β, device=device)


"""
    Free energy: F = -1/(βL) lnZ
"""
function free_energy(ψl::CMPS, ψr::CMPS, W::CMPO, β::Real; device=:cpu)
    res = log_overlap(ψl, W * ψr, β, device=device) - log_overlap(ψl, ψr, β, device=device)
    return -res/β
end
free_energy(ψ::CMPS, W::CMPO, β::Real; device=:cpu) = free_energy(ψ, ψ, W, β, device=device)


"""
    Energy density: E = -∂lnZ/∂β
"""
function energy(ψl::CMPS, ψr::CMPS, W::CMPO, β::Real; device=:cpu)
    K = ψl * W * ψr 
    H = ψl * ψr 
    res = thermal_average(K, ψl, ψr, W, β, device=device) - thermal_average(H, ψl, ψr, β, device=device)
    return res
end
energy(ψ::CMPS, W::CMPO, β::Real; device=:cpu) = energy(ψ, ψ, W, β, device=device)


"""
    Specific heat: Cv = -β^2 ∂E/∂β
"""
function specific_heat(ψl::CMPS, ψr::CMPS, W::CMPO, β::Real; 
                        method::Symbol = :adiff,
                        device = :cpu)
    if method == :adiff
        K = ψl * W * ψr 
        H = ψl * ψr 
        K2 = K * K
        H2 = H * H
        c = thermal_average(K2, ψl, ψr, W, β, device=device) - thermal_average(K, ψl, ψr, W, β, device=device)^2
        c -= thermal_average(H2, ψl, ψr, β, device=device) - thermal_average(H, ψl, ψr, β, device=device)^2
    elseif method == :ndiff
        minus_e = b -> -energy(ψl, ψr, W, b, device=device)
        c = central_fdm(5, 1)(minus_e, β)
    else @error "method should be :adiff or :ndiff"
    end
    return β^2 * c
end
specific_heat(ψ::CMPS, W::CMPO, β::Real; 
                method::Symbol = :adiff,
                device=:cpu) = specific_heat(ψ, ψ, W, β, method = method, device=device)


"""
    Entropy: S = β×(E - F)
"""
function entropy(ψl::CMPS, ψr::CMPS, W::CMPO, β::Real; device=:cpu)
    s = energy(ψl, ψr, W, β, device=device) - free_energy(ψl, ψr, W, β, device=device)
    return β*s
end
entropy(ψ::CMPS, W::CMPO, β::Real; device=:cpu) = entropy(ψ, ψ, W, β, device=device)

#end  # module ThermaldynamicQuanties
