#module PhysicalModels
#include("Setup.jl")

"""
NN Transvers field Ising model
    H = ∑ J Zi Zj + ∑ Γ Xi
"""
function TFIsing(J::Real, Γ::Real; field::Symbol=:N, η::Float64 = 1.e-2)
    if field == :N
        h = zeros(2,2)
    else
        h = η .* pauli(field)
    end
    return CMPO(Γ*pauli(:x)+h, √J*pauli(:z), √J*pauli(:z), zeros(2,2))
end

"""
XY model
    H = ∑ Xi Xj + Yi Yj
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ (S+ S+ + S-S-)
"""
function XYmodel()
    sp = pauli(:+); sm = pauli(:-);
    L = zeros(2,2,2)
    L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm
    R = zeros(2,2,2)
    R[:,:,1] = 1/√2 * sp ; R[:,:,2] = 1/√2 * sm
    Q = zeros(2,2)
    P = zeros(2,2,2,2)
    return CMPO(Q,R,L,P)
end


"""
Heisenberg XXZ model
    H = ∑ Xi Xj + Yi Yj + Δ Zi Zj
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ (S+ S+ + S-S-) - ∑ Δ Zi Zj
"""
function XXZmodel(Δ::Real)
    if Δ == 0
        return XYmodel()
    else
        sp = pauli(:+); sm = pauli(:-); sz = 0.5 * pauli(:z)
        L = zeros(2,2,3)
        L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm; L[:,:,3] = √Δ * sz
        R = zeros(2,2,3)
        R[:,:,1] = 1/√2 * sp ; R[:,:,2] = 1/√2 * sm; R[:,:,3] = √Δ * sz
        Q = zeros(2,2)
        P = zeros(2,2,3,3)
        return CMPO(Q,R,L,P)
    end
end

"""
Heisenberg model
    H = J ∑ Xi Xj + Yi Yj + Zi Zj
AFM case: J = 1.0, after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ (S+ S+ + S-S-) - ∑ Zi Zj
"""
function HeisenbergModel(; J::Real=1.0)
    if J > 0
        return XXZmodel(1.0)
    else
        @error "Not support yet!"
    end
end







#end module PhysicalModels
