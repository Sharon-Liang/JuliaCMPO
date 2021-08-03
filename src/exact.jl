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

function free_energy(J::Real, Γ::Real, β::Real; L = 10000)
    f = 0.
    for n = 1:L
        k = 2*π/L * n
        ϵ = energy_density(k, J, Γ)
        f += logaddexp(β*ϵ/2, -β*ϵ/2)
    end
    return -f/(β*L)
end

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

function critical_zz_cor(τ::Real, β::Real;J::Real=1.0)
    g0 = 0.858714569; zc = J^(-1/4)
    T = 1/β
    fac = zc * 2^(-1/4)* T^(1/4) * g0
    down = sin(π*T*τ)
    down = down^(1/4)
    return fac/down |> real
end

function critical_zz_sus(n::Integer, β::Real; J::Real=1.0)
    g0 = 0.858714569; zc = J^(-1/4)
    T = 1/β
    fac = 4/3 * zc * g0 * β^(3/4)
    ωn = 2π * n/β  #bosion
    res = exp(1im*n*π) / beta(7/8 + n,7/8 - n)
    return fac*res |> real
end


function critical_zz_chi(ω::Real, β::Real; J::Real=1.0)
    # Ref: PRX.4.031008(2014) A6
    g0 = 0.858714569; zc = J^(-1/4)
    T = 1/β
    up = zc * g0 * β^(3/4)
    down = 2^(1/4) * √π * gamma(1/8) * gamma(5/8)
    fac = up/down
    res = sinh(ω/(2*T)) * abs(gamma(1/8 - 1im*ω/(2π*T)))^2
    return fac*res
end

function structure_factor(J::Real, Γ::Real, β::Real)
    Δ = 2*√(J^2 - Γ^2)
    if Δ==0
        g0 = 0.858714569; zc = J^(-1/4)
        up = zc * g0 * β^(3/4) * gamma(1/8)
        down = 2^(1/4) * √π * gamma(5/8)
        res = up/down
    else
        error("Not support yet.")
    end
end
#end
