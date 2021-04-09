# exact solutions of TFIsing model
function energy_density(k::Real, J::Real, Γ::Real)
    eng = 2*sqrt(J^2+ Γ^2 - 2*J*Γ*cos(k))
    return eng
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

function critical_zz_cor(β::Real, τ::Real, x::Real ;J::Real=1.0)
    g0 = 0.858714569; zc = (2*J)^(-1/4); c = 2*J
    T = 1/β
    fac = zc * T^(1/4) * g0
    down = sin(π*T*(τ-1im*x/c)) * sin(π*T*(τ+1im*x/c))
    down = down^(1/8)
    return fac/down |> real
end
