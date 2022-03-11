#module PhysicalModels
#include("Setup.jl")

"""
    Ising type CMPO: H = -J ôl ôr block
"""
function Ising_CMPO(J::Real, ol::AbstractArray, or::AbstractArray, wid::Integer = 1)
    sgn = sign(J); val =√(abs(J))
    o2 = zeros(2, 2)
    i2 = Matrix{Float64}(I, 2, 2)
    #construct L and R
    ol = val * ol; or = sgn * val * or
    L = ol; R = or
    for i = 2: wid
        R = cat(R, o2, dims = 3)
        i == wid ? L = cat(L, ol, dims=3) : L = cat(L, o2, dims = 3)
    end

    #construct P 
    P = []
    for i = 1: wid
        pcol = o2
        for j = 2 :wid
            j == i+1 ? pcol = cat(pcol, i2, dims = 3) : pcol = cat(pcol, o2, dims = 3)
        end
        i == 1 ? P = pcol : P = cat(P, pcol, dims = 4)
    end
        
    return CMPO(zeros(2,2), R, L, P)
end


"""
    General form of W matrix: W = [0 I; I 0]
"""
function generalUt(vir_dim::Integer)
    zero_block = zeros(vir_dim, vir_dim)
    eye_block  = Matrix{Float64}(I, vir_dim, vir_dim)
    col1 = cat(zero_block, eye_block, dims = 1)
    col2 = cat(eye_block, zero_block, dims = 1)
    return cat(col1, col2, dims = 2), 2*vir_dim
end


"""
    Expand cmpo
"""
function expand(o::CMPO)
    R = √(0.5) * cat(o.R, o.L, dims = 3)
    L = √(0.5) * cat(o.L, o.R, dims = 3)

    if length(size(o.P)) == 2
        oP = reshape(o.P, size(o.P)[1], size(o.P)[2], 1, 1)
    else
        oP = o.P
    end

    pl = cat(oP, zeros(size(oP)), dims = 3)
    pr = cat(zeros(size(oP)), ein"ijkl->ijlk"(oP), dims = 3)
    P = cat(pl, pr, dims = 4)
    
    return CMPO(o.Q, R, L, P)
end


"""
    cat CMPO blocks
"""
function cat(o1::CMPO, o2::CMPO)
    pd = size(o1.Q)[1]
    length(size(o1.P)) == 2 ? vd1 = 1 : vd1 = size(o1.P)[3]
    length(size(o2.P)) == 2 ? vd2 = 1 : vd2 = size(o2.P)[3]

    Q = zeros(pd, pd)
    R = cat(o1.R, o2.R, dims = 3)
    L = cat(o1.L, o2.L, dims = 3)

    pl = cat(o1.P, zeros(pd, pd, vd2, vd1), dims = 3)
    pr = cat(zeros(pd, pd, vd1, vd2), o2.P, dims = 3)
    P = cat(pl, pr, dims = 4)

    return CMPO(Q, R, L, P)
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
    Tmatrix = CMPO(Γ*pauli(:x)+h, √J*pauli(:z), √J*pauli(:z), zeros(2,2))
    Ut = Matrix{Float64}(I, 1, 1)
    return PhysModel(Tmatrix, 2, 1, Ut)
end


"""
NN XY model
    H = ∑ Xi Xj + Yi Yj
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ (S+ S+ + S-S-)
"""
function XYmodel()
    sp = pauli(:+); sm = pauli(:-);
    L = zeros(2, 2, 2)
    L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm
    R = zeros(2, 2, 2)
    R[:,:,1] = 1/√2 * sp ; R[:,:,2] = 1/√2 * sm
    Q = zeros(2, 2)
    P = zeros(2, 2, 2, 2)
    Tmatrix = CMPO(Q,R,L,P)
    Ut = Matrix(1.0I, 2, 2)
    return PhysModel(Tmatrix, 2, 2, Ut)
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
        L = zeros(2, 2, 3)
        L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm; L[:,:,3] = √Δ * sz
        R = zeros(2, 2, 3)
        R[:,:,1] = 1/√2 * sp ; R[:,:,2] = 1/√2 * sm; R[:,:,3] = √Δ * sz
        Q = zeros(2, 2)
        P = zeros(2, 2, 3, 3)
        Tmatrix = CMPO(Q,R,L,P)
        Ut = Matrix(1.0I, 3, 3) 
        return PhysModel(Tmatrix, 2, 3, Ut)
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
function TFIsing_2D_helical(J::Real, Γ::Real, wid::Integer = 1; field::Symbol=:N, η::Float64 = 1.e-2)
    if field == :N
        h = zeros(2,2)
    else
        h = η .* pauli(field)
    end
    Q = Γ*pauli(:x)+h
    T = Ising_CMPO(J, pauli(:z), pauli(:z), wid)
    Tmatrix = CMPO(Q, T.R, T.L, T.P) |> expand
    Ut, vir_dim = generalUt(wid)
    return PhysModel(Tmatrix, 2, vir_dim, Ut)
end


"""
2D XY model, helical boundary condition
    H = ∑ [Xi X_(i+1) +Xi X_(i+W) + Yi Y_(i+1) +Yi Y_(i+W)]
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ [S+_i S+_(i+1) + S-_iS-_(i+1) + S+_i S+_(i+W) + S-_iS-_(i+W)]
"""
function XYmodel_2D_helical(wid::Integer = 1)
    phy_dim = 2; vir_dim = 2*wid
    sp = pauli(:+); sm = pauli(:-);
    Tp = Ising_CMPO(0.5, sp, sp, wid)
    Tm = Ising_CMPO(0.5, sm, sm, wid)
    Tmatrix = cat(Tp, Tm) |> expand
    Ut, vir_dim = generalUt(vir_dim)
    return PhysModel(Tmatrix, phy_dim, vir_dim, Ut)
end


"""
2D Heisenberg XXZ model helical boundary condition
    H = ∑ [Xi X_(i+1) +Xi X_(i+W) + Yi Y_(i+1) +Yi Y_(i+W)] + Δ [Zi Z_(i+1) +Zi Z_(i+W)]
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ [S+_i S+_(i+1) + S-_iS-_(i+1) + S+_i S+_(i+W) + S-_iS-_(i+W)] - Δ [Zi Z_(i+1) +Zi Z_(i+W)]
"""
function XXZmodel_2D_helical(Δ::Real, wid::Integer = 1)
    phy_dim = 2; vir_dim = 3*wid
    sp = pauli(:+); sm = pauli(:-); sz = 0.5 * pauli(:z)
    Tp = Ising_CMPO(0.5, sp, sp, wid)
    Tm = Ising_CMPO(0.5, sm, sm, wid)
    Tz = Ising_CMPO(Δ, sz, sz, wid)
    Tmatrix = cat(cat(Tp, Tm), Tz) |> expand
    Ut, vir_dim = generalUt(vir_dim)
    return PhysModel(Tmatrix, phy_dim, vir_dim, Ut)
end

#end module PhysicalModels
