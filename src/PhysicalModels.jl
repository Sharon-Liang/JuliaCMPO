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
function XYmodel(; J::Real = 1.0)
    sgn = sign(J) 
    sp = pauli(:+); sm = pauli(:-);
    L = zeros(2,2,2)
    L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm
    R = zeros(2,2,2)
    R[:,:,1] = -sgn*1/√2 * sm ; R[:,:,2] = -sgn*1/√2 * sp
    Q = zeros(2,2)
    P = zeros(2,2,2,2)
    return cmpo(Q,R,L,P)
end
"""
function HeisenbergModel(; J::Real=1.0)
    #AFM Heisenberg model
    sp = pauli(:+); sm = pauli(:-);
    LR = zeros(2,2,3)
    LR[:,:,1] = 1/√2 * sp ; LR[:,:,2] = 1/√2 * sm; LR[:,:,3] = pauli(:z)
    Q = zeros(2,2)
    P = zeros(2,2,3,3) 
    return cmpo(Q,LR,LR,P)
end

"""
function XXZmodel(Jx::Real, Jz::Real)
    if Jx == Jz 
        return HeisenbergModel(J=Jz)
    elseif Jz == 0
        return XYmodel(J=Jx)
    elseif Jx == 0
        return TFIsing(Jz, 0.0)
    else
        sgnz = sign(Jz) ; sgnx = sign(Jx)
        Jz = abs(Jz); Jx = abs(Jx)
        sp = pauli(:+); sm = pauli(:-);
        L = zeros(2,2,3)
        L[:,:,1] = √(Jx/2) * sp ; L[:,:,2] = √(Jx/2) * sm; L[:,:,3] = √(Jz) *pauli(:z)
        R = zeros(2,2,3)
        R[:,:,1] = -sgn * √(Jx/2) * sm ; R[:,:,2] = -sgn * √(Jx/2) * sp; R[:,:,3] = -sgn * √(Jz) *pauli(:z)
        Q = zeros(2,2)
        P = zeros(2,2,3,3)
        return cmpo(Q,R,L,P)
    end
end
"""

#end module PhysicalModels
