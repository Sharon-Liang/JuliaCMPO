#module PhysicalModels
#include("Setup.jl")
struct PhysModel
    Tm::CMPO # local Transfer matrix
    pD::Int64 #pD: bond dimension of physical legs
    vD::Int64 # vD+1: virtual bond dimension of the horizontal legs
end

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
    return PhysModel(CMPO(Γ*pauli(:x)+h, √J*pauli(:z), √J*pauli(:z), zeros(2,2)),2,1)
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
    return PhysModel(CMPO(Q,R,L,P),2,2)
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
        return PhysModel(CMPO(Q,R,L,P),2,3)
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


"""
2D NN Transvers field Ising model,helical boundary condition
    H = ∑ J [Zi Z_(i+1) +Zi Z_(i+W)]  + ∑ Γ Xi
"""
function TFIsing(J::Real, Γ::Real, W::Int64; field::Symbol=:N, η::Float64 = 1.e-2)
    if W == 1
        return TFIsing(J,Γ)
    else
        if field == :N
            h = zeros(2,2)
        else
            h = η .* pauli(field)
        end
        Q = Γ*pauli(:x)+h
        L = zeros(2,2,W)
        L[:,:,1] = √J*pauli(:z) ; L[:,:,W] = √J*pauli(:z)
        R = zeros(2,2,W)
        R[:,:,1] = √J*pauli(:z)
        P = zeros(2,2,W,W)
        for i=2:W
            P[:,:,i,i-1] = Matrix{Float64}(I,2,2)
        end
        return PhysModel(CMPO(Q,R,L,P),2,W)
    end
end


"""
2D XY model, ,helical boundary condition
    H = ∑ [Xi X_(i+1) +Xi X_(i+W) + Yi Y_(i+1) +Yi Y_(i+W)]
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ [S+_i S+_(i+1) + S-_iS-_(i+1) + S+_i S+_(i+W) + S-_iS-_(i+W)]
"""
function XYmodel(W::Int64)
    if W == 1
        return XYmodel()
    else
        sp = pauli(:+); sm = pauli(:-);
        L = zeros(2,2,2W)
        L[:,:,1] = 1/√2 * sp ; L[:,:,W] = 1/√2 * sp
        L[:,:,W+1] = 1/√2 * sm; L[:,:,2W] = 1/√2 * sm
        R = zeros(2,2,2W)
        R[:,:,1] = 1/√2 * sp ; R[:,:,W+1] = 1/√2 * sm
        Q = zeros(2,2)
        P = zeros(2,2,2W,2W)
        for j=0:1, i=2:W
            P[:,:,j*W+i,j*W+i-1] = Matrix{Float64}(I,2,2)
        end
        return PhysModel(CMPO(Q,R,L,P),2,2W)
    end
end

"""
2D Heisenberg XXZ model helical boundary condition
    H = ∑ [Xi X_(i+1) +Xi X_(i+W) + Yi Y_(i+1) +Yi Y_(i+W)] + Δ [Zi Z_(i+1) +Zi Z_(i+W)]
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ [S+_i S+_(i+1) + S-_iS-_(i+1) + S+_i S+_(i+W) + S-_iS-_(i+W)] - Δ [Zi Z_(i+1) +Zi Z_(i+W)]
"""
function XXZmodel(Δ::Real, W::Int)
    if W == 1
        return XXZmodel(Δ)
    else
        if Δ == 0
            return XYmodel(W)
        else
            sp = pauli(:+); sm = pauli(:-); sz = 0.5 * pauli(:z)
            L = zeros(2,2,3W)
            L[:,:,1] = 1/√2 * sp ; L[:,:,W] = 1/√2 * sp
            L[:,:,W+1] = 1/√2 * sm; L[:,:,2W] = 1/√2 * sm
            L[:,:,2W+1] = √Δ * sz; L[:,:,3W] = √Δ * sz
            R = zeros(2,2,3W)
            R[:,:,1] = 1/√2 * sp 
            R[:,:,W+1] = 1/√2 * sm
            R[:,:,2W+1] = √Δ * sz
            Q = zeros(2,2)
            P = zeros(2,2,3W,3W)
            for j=0:2, i=2:W
                P[:,:,j*W+i,j*W+i-1] = Matrix{Float64}(I,2,2)
            end
            return PhysModel(CMPO(Q,R,L,P),2,3W)
        end
    end
end


#end module PhysicalModels
