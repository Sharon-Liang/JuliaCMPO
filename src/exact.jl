#using LinearAlgebra
#using StatsFuns, SpecialFunctions, HCubature
#import cMPO:free_energy, energy, specific_heat, entropy

# exact solutions of TFIsing model
function energy_density(k::Real, J::Real, Γ::Real)
    eng = 2*sqrt(J^2+ Γ^2 - 2*J*Γ*cos(k))
    return eng
end

function energy(J::Real, Γ::Real, β::Real; L = 10000)
    e = 0.
    for n = 1:L
        k = 2*π/L * n
        ϵ = energy_density(k, J, Γ)
        e += ϵ *(logistic(-β*ϵ) - 0.5)
    end
    return e/L
end

function fk(k::Real,J::Real, Γ::Real, β::Real)
    ϵ = energy_density(k, J, Γ)
    return -logaddexp(β*ϵ/2, -β*ϵ/2)/β
end

function free_energy(J::Real, Γ::Real, β::Real; err::Float64=eps())
    res = hquadrature(k->fk(k,J,Γ,β)/2π, 0, 2π,rtol=err)
    return res
end

"""
function free_energy(J::Real, Γ::Real, β::Real; L::Integer = 10000)
    f = 0.
    for n = 1:L
        k = 2*π/L * n
        ϵ = energy_density(k, J, Γ)
        f += logaddexp(β*ϵ/2, -β*ϵ/2)
    end
    return -f/(β*L)
end
"""

function entropy(J::Real, Γ::Real, β::Real; L = 10000)
    s = energy(J, Γ, β; L) - free_energy(J, Γ, β; L)
    return β*s
end

function specific_heat(J::Real, Γ::Real, β::Real; L = 10000)
    c = 0.
    for n = 1:L
        k = 2*π/L * n
        ϵ = energy_density(k, J, Γ)
        c += ϵ^2*(logistic(-β*ϵ) - logistic(-β*ϵ)^2)                                                                                                                 
    end
    return β^2 * c/L
end

function ave_sx(J::Real, Γ::Real, β::Real; L = 100)
    num = 0.
    L2 = 2*L
    for n = -L+1:L
        k = 2*π/L2 * n
        θk = atan(-J*sin(k), -J*cos(k)+Γ) # J and Γ with signs
        uk = cos(θk/2) ; vk = sin(θk/2)
        ϵ = energy_density(k, J, Γ)
        num += uk^2 * logistic(-β*ϵ) + vk^2 * logistic(β*ϵ)
    end
    return  1 - 2*num/L2
end

function correlation_2time(J::Real,τ::Real, β::Real;Γ::Real=1.0)
    if J/Γ == 1.
        g0 = 0.858714569; zc = J^(-1/4)
        T = 1/β
        num = zc * 2^(-1/4)* T^(1/4) * g0
        den = sin(π*T*τ)
        den = den^(1/4)
        return num/den |> real
    elseif J/Γ == 0.
        return cosh(-β*Γ + 2Γ*τ)/cosh(-β*Γ)
    else
        @error "No exact results. J/Γ should be 1. or 0."
    end
end

function spectral_density(J::Real, ω::Real, β::Real; η::Float64=0.05,Γ::Real=1.)
    # 2Imχ(ω)
    if J/Γ == 0.
        return 2π*sinh(Γ*β)/cosh(Γ*β)*(delta(ω-2Γ,η) - delta(ω+2Γ,η))
    elseif J/Γ == 1.
        # Ref: PRX.4.031008(2014) A6
        g0 = 0.858714569; zc = J^(-1/4)
        T = 1/β
        up = zc * g0 * β^(3/4)
        down = 2^(1/4) * √π * gamma(1/8) * gamma(5/8)
        fac = up/down
        res = sinh(ω/(2*T)) * abs(gamma(1/8 - 1im*ω/(2π*T)))^2
        return 2*fac*res
    else
        @error "No exact results. J/Γ should be 1. or 0."
    end
end

function susceptibility(J::Real, ω::Real, β::Real; η::Float64=0.05,Γ::Real=1.)
    return spectral_density(J, ω, β; η=η,Γ=Γ)
end

function structure_factor(J::Real, ω::Real, β::Real; η::Float64=0.05,Γ::Real=1.)
    if J/Γ == 0.
        return 1/(1-exp(-β*ω))*spectral_density(J,ω,β,Γ=Γ,η=η)
    elseif J/Γ == 1.
        if ω == 0.
            g0 = 0.858714569; zc = J^(-1/4)
            up = zc * g0 * β^(3/4) * gamma(1/8)
            down = 2^(1/4) * √π * gamma(5/8)
            return up/down
        else
            @error "only ω=0 result avaliable"
        end
    else
        @error "No exact results. J/Γ should be 1. or 0."
    end
end

#end
