""" Exact solutions """

"""NN Transverse Field Ising Model """
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
        θk = atan(J*sin(k), J*cos(k)-Γ)
        uk = cos(θk/2) ; vk = sin(θk/2)
        ϵ = energy_density(k, J, Γ)
        num += uk^2 * logistic(-β*ϵ) + vk^2 * logistic(β*ϵ)
    end
    return  1 - 2*num/L2
end
