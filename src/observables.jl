#=
### *Thermaldynamic quanties*
=#

""" 
    free_energy(ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)

Calculate the total free_energy:
```math
F = -\\frac{1}{β}lnZ
```
where ``β = 1/T`` is the inverse temperature.
"""
function free_energy(ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)
    res = log_overlap(ψl, Tₘ * ψr, β) - log_overlap(ψl, ψr, β)
    return -res/β
end

free_energy(ψ::CMPS, Tₘ::CMPO, β::Real) = free_energy(ψ, ψ, Tₘ, β)



"""
    thermal_average(Ô::AbstractMatrix, ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)

The thermal average of local opeartors ``Ô`` with respect to ``K = ψl * Tₘ * ψr``,
```math
    ⟨Ô⟩ = \frac{Tr(Ô e^{-βK})}{Z}
```

Note: The input operator ``Ô`` could be with the same physical bond dimension as the cMPO ``Tₘ`` or with the same size as ``K``.
"""
function thermal_average(Ô::AbstractMatrix, ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)
    K = ψl * Tₘ * ψr |> symmetrize
    if size(Ô) == size(Tₘ.Q)
        eye = oneunit(ψl.Q)
        Ô = eye ⊗ Ô ⊗ eye
    end

    @assert size(Ô) == size(K)
    
    e, v = eigensolver(K)
    e₀ = minimum(e)
    e = e .- e₀
    Λ = @. exp(-β * e)
    ave = diag(v' * Ô * v)

    den = sum(Λ)
    num = reduce(+, map(*, Λ, ave))
    return num/den   
end

thermal_average(Ô::AbstractArray, ψ::CMPS, Tₘ::CMPO, β::Real) = thermal_average(Ô, ψ, ψ, Tₘ, β)



"""
    thermal_average(Ô::AbstractArray, ψl::CMPS, ψr::CMPS, β::Real)

The thermal average of local opeartors ``Ô`` with respect to ``K = ψl * ψr``,
```math
    ⟨Ô⟩ = \\frac{Tr(Ô e^{-βK})}{Z}
```
Note: The input operator ``Ô`` could be with the same size as ``K``.
"""
function thermal_average(Ô::AbstractArray, ψl::CMPS, ψr::CMPS, β::Real)
    K = ψl * ψr |> symmetrize
    @assert size(Ô) == size(K)
    
    e, v = eigensolver(K)
    e₀ = minimum(e) 
    e = e .- e₀
    Λ = @. exp(-β * e)
    ave = diag(v' * Ô * v)

    den = sum(Λ)
    num = reduce(+, map(*, Λ, ave))
    return num/den   
end

thermal_average(Ô::AbstractArray, ψ::CMPS, β::Real) = thermal_average(Ô, ψ, ψ, β)



"""
    energy(ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)

Calculate the energy density: ``E = -∂lnZ/∂β``.
"""
function energy(ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)
    K = ψl * Tₘ * ψr 
    H = ψl * ψr 
    res = thermal_average(K, ψl, ψr, Tₘ, β) - thermal_average(H, ψl, ψr, β)
    return res
end

energy(ψ::CMPS, Tₘ::CMPO, β::Real) = energy(ψ, ψ, Tₘ, β)


"""
    entropy(ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)

Calculate the entropy: ``S = β × (E - F)``.
"""
function entropy(ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)
    E = energy(ψl, ψr, Tₘ, β)
    F = free_energy(ψl, ψr, Tₘ, β)
    return β*(E-F)
end
entropy(ψ::CMPS, Tₘ::CMPO, β::Real) = entropy(ψ, ψ, Tₘ, β)



"""
    specific_heat(ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)

Calculate specific heat: ``Cᵥ = -β² ∂E/∂β``
"""
function specific_heat(ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)
    K = ψl * Tₘ * ψr 
    H = ψl * ψr 

    K² = K * K
    H² = H * H

   res = thermal_average(K², ψl, ψr, Tₘ, β) - thermal_average(K, ψl, ψr, Tₘ, β)^2
   res -= thermal_average(H², ψl, ψr, β) - thermal_average(H, ψl, ψr, β)^2

   return β^2 * res
end
specific_heat(ψ::CMPS, Tₘ::CMPO, β::Real) = specific_heat(ψ, ψ, Tₘ, β)




#=
### *Correlations*
=#
"""
    correlation_2time(τ::Number, Ô₁::AbstractMatrix,Ô₂::AbstractMatrix,ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)

Calculate the local two-time correlation functions: ``C(τ) = -G(τ) = ⟨A(τ)B(0)⟩.``
"""
function correlation_2time(τ::Number, Ô₁::AbstractMatrix,Ô₂::AbstractMatrix, ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real)
    K = ψl * Tₘ * ψr  |> symmetrize
    e, v = eigensolver(K)

    e₀ = minimum(e)
    e = e .- e₀

    if size(Ô) == size(Tₘ.Q)
        eye = oneunit(ψl.Q)
        Ô₁ = eye ⊗ Ô₁ ⊗ eye
        Ô₂ = eye ⊗ Ô₂ ⊗ eye
    end
    @assert size(Ô₁) == size(K)
    @assert size(Ô₂) == size(K)

    Ô₁ = v' * Ô₁ * v
    Ô₂ = v' * Ô₂ * v

    den = exp.(-β * e) |> sum
    num = 0.0
    for i in eachindex(e), j in eachindex(e)
        num += exp(-β*e[i] + τ*(e[i] - e[j])) * Ô₁[i,j] * Ô₂[j,i]
    end
    return num/den
end

correlation_2time(τ::Number, Ô₁::AbstractMatrix,Ô₂::AbstractMatrix, ψ::CMPS, Tₘ::CMPO, β::Real) = correlation_2time(τ, Ô₁, Ô₂, ψ, ψ, Tₘ, β) 


"""
    _kernel_Lehmann(z::Number, Eₙ::Real, Eₘ::Real, β:: Real, type::OperatorType = Bose)

Calculate the kernal function in Lehmann representation:
```math
    K(z, Eₙ, Eₘ, β, λ) = \\frac{exp(-βEₙ) - λ exp(-βEₘ)}{z - Eₘ + Eₙ}
```
for Masubara frequency Green's function, where ``λ = 1`` when ``T=Bose`` and  ``λ = -1`` when ``T=Fermi``.
"""
function _kernel_Lehmann(z::Number, Eₙ::Real, Eₘ::Real, β::Real, type::OperatorType = Bose)
    type == Bose ? λ = 1.0 : λ = -1.0
    if type == Bose && z == 0 && abs(Eₙ - Eₘ) < 1.e-10
        return -β * exp(-β * Eₘ)
    else
        num = exp(-β*Eₙ) - λ * exp(-β*Eₘ)
        den = z - Eₘ + Eₙ
        return num/den
    end
end



"""
    retarded_GF(z::Number, Ô₁::AbstractMatrix,Ô₂::AbstractMatrix, ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real, type::OperatorType = Bose)

Calculate the Lehmann representation of retarded Greens function `G(z)`.
"""
function retarded_GF(z::Number, Ô₁::AbstractMatrix,Ô₂::AbstractMatrix, ψl::CMPS, ψr::CMPS, Tₘ::CMPO, β::Real, type::OperatorType = Bose)

    K = ψl * Tₘ * ψr  |> symmetrize
    e, v = eigensolver(K)

    e₀ = minimum(e)
    e = e .- e₀

    if size(Ô) == size(Tₘ.Q)
        eye = oneunit(ψl.Q)
        Ô₁ = eye ⊗ Ô₁ ⊗ eye
        Ô₂ = eye ⊗ Ô₂ ⊗ eye
    end
    @assert size(Ô₁) == size(K)
    @assert size(Ô₂) == size(K)

    Ô₁ = v' * Ô₁ * v
    Ô₂ = v' * Ô₂ * v

    den = exp.(-β * e) |> sum
    num = 0.0

    for i in eachindex(e), j in eachindex(e)
        num += Ô₁[i,j] * Ô₂[j,i] * _kernel_Lehmann(z, e[i], e[j], β, type)
    end
    return num/den
end


retarded_GF(z::Number, Ô₁::AbstractMatrix,Ô₂::AbstractMatrix, ψ::CMPS, Tₘ::CMPO, β::Real, type::OperatorType = Bose) = retarded_GF(z, Ô₁,Ô₂, ψ, ψ, Tₘ, β, type)




