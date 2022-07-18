#module PhysicalModels
"""
    Data structure of a physical model 
    PhysModel:
        `Tmatrix`: local transfer matrix
        `phy_dim`: bond dimension of physical legs
        `vir_dim`: virtual bond dimension of the horizontal legs
        `Ut`: unitary transformation that makes: `Ut^† Tmatrix Ut = transpose(Tmatrix)`
              Note that in general, `Ut` not necessarily be unitary, but be invertible, 
              however, in the limiting cases we know right now, they are unitary.
"""
@with_kw struct PhysModel{Tm<:AbstractCMPO, 
                          Tu<:Union{AbstractMatrix, Nothing}}
    Tmatrix::Tm 
    phy_dim::Int64  
    vir_dim::Int64  
    Ut::Tu
end


"""
    Ising type CMPO: H = -J ôl ôr block
"""
function Ising_CMPO(J::Real, ol::AbstractArray{T}, or::AbstractArray{T}, wid::Integer = 1) where {T}
    sgn = sign(J); val =√(abs(J))
    o2 = zeros(2, 2)
    i2 = Matrix{T}(I, 2, 2)
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
        
    return CMPO(zeros(T,2,2), R, L, P)
end


"""
    General form of `Ut` matrix: `Ut = [0 I; I 0]`.
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
function expand_cmpo(o::AbstractCMPO{T,S,U,V}) where {T,S,U,V}
    R = √(0.5) * cat(o.R, o.L, dims = 3)
    L = √(0.5) * cat(o.L, o.R, dims = 3)

    if V <: AbstractMatrix
        oP = reshape(o.P, size(o.P,1), size(o.P,2), 1, 1)
    else
        oP = o.P
    end

    pl = cat(oP, zeros(size(oP)), dims = 3)
    pr = cat(zeros(T,size(oP)), ein"ijkl->ijlk"(oP), dims = 3)
    P = cat(pl, pr, dims = 4)
    
    return CMPO(o.Q, R, L, P)
end


"""
    cat CMPO blocks
"""
function cat(o1::AbstractCMPO{T,S1,U1,V1}, 
             o2::AbstractCMPO{T,S2,U2,V2}
             ) where {T,S1,U1,V1, S2,U2,V2}
    pd = size(o1.Q,1)
    V1<:AbstractMatrix ? vd1 = 1 : vd1 = size(o1.P,3)
    V2<:AbstractMatrix ? vd2 = 1 : vd2 = size(o2.P,3)

    Q = zeros(T, pd, pd)
    R = cat(o1.R, o2.R, dims = 3)
    L = cat(o1.L, o2.L, dims = 3)

    pl = cat(o1.P, zeros(T, pd, pd, vd2, vd1), dims = 3)
    pr = cat(zeros(T, pd, pd, vd1, vd2), o2.P, dims = 3)
    P = cat(pl, pr, dims = 4)

    return CMPO(Q, R, L, P)
end


"""
NN Transvers field Ising model
    H = -∑ J Zi Zj - ∑ Γ Xi
"""
function TFIsing(J::Real, Γ::Real; 
            testfield::Union{PauliMatrixName,Nothing} = nothing, 
            η::Float64 = 1.e-2)
    testfield === nothing ? h = zeros(2,2) : h = η .* pauli(testfield)
    Tmatrix = CMPO(Γ*pauli(PX)+h, √J*pauli(PZ), √J*pauli(PZ), zeros(2,2))
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
    sp = pauli(PPlus); sm = pauli(PMinus);
    L = zeros(2, 2, 2)
    L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm
    R = zeros(2, 2, 2)
    R[:,:,1] = 1/√2 * sp ; R[:,:,2] = 1/√2 * sm
    Q = zeros(2, 2)
    P = zeros(2, 2, 2, 2)
    Tmatrix = CMPO(Q,R,L,P)
    Ut = Matrix{Float64}(I, 2, 2)
    return PhysModel(Tmatrix, 2, 2, Ut)
end


"""
Heisenberg XXZ model
    H = ∑ Xi Xj + Yi Yj + Δ Zi Zj
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ (S+ S+ + S-S-) - ∑ Δ Zi Zj
"""

"""
function XXZmodel_tw(Δ::Real)
    if Δ == 0
        return XYmodel()
    else
        sp = pauli(PPlus); sm = pauli(PMinus); sz = 0.5 * pauli(PZ)
        L = zeros(2, 2, 3)
        L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm; L[:,:,3] = √Δ * sz
        R = zeros(2, 2, 3)
        R[:,:,1] = 1/√2 * sm ; R[:,:,2] = 1/√2 * sp; R[:,:,3] = -√Δ * sz
        Q = zeros(2, 2)
        P = zeros(2, 2, 3, 3)
        Tmatrix = CMPO(Q,R,L,P)
        Ut = [0 1 0; 1 0 0; 0 0 -1] 
        return PhysModel(Tmatrix, 2, 3, Ut)
    end
end
"""

function XXZmodel(Δ::Real)
    if Δ == 0
        return XYmodel()
    else
        sp = pauli(PPlus); sm = pauli(PMinus); sz = 0.5 * pauli(PZ)
        L = zeros(2, 2, 3)
        L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm; L[:,:,3] = √Δ * sz
        R = zeros(2, 2, 3)
        R[:,:,1] = 1/√2 * sp ; R[:,:,2] = 1/√2 * sm; R[:,:,3] = √Δ * sz
        Q = zeros(2, 2)
        P = zeros(2, 2, 3, 3)
        Tmatrix = CMPO(Q,R,L,P)
        Ut = Matrix{Float64}(I, 3, 3) 
        return PhysModel(Tmatrix, 2, 3, Ut)
    end
end


"""
2D NN Transvers field Ising model,helical boundary condition
    H = -∑ J [Zi Z_(i+1) +Zi Z_(i+W)]  - ∑ Γ Xi
"""
function TFIsing_2D_helical(J::Real, Γ::Real, wid::Integer = 1; 
                            testfield::Union{PauliMatrixName,Nothing} = nothing, 
                            η::Float64 = 1.e-2, 
                            expand::Bool = false)
    testfield === nothing ? h = zeros(2,2) : h = η .* pauli(field)
    phy_dim = 2; vir_dim = wid
    Q = Γ*pauli(PX)+h
    T = Ising_CMPO(J, pauli(PZ), pauli(PZ), wid)
    Tmatrix = CMPO(Q, T.R, T.L, T.P)
    Ut = nothing
    if expand
        Tmatrix = Tmatrix |> expand_cmpo
        Ut, vir_dim = generalUt(vir_dim)
    end
    return PhysModel(Tmatrix, phy_dim, vir_dim, Ut)
end


"""
2D XY model, helical boundary condition
    H = ∑ [Xi X_(i+1) +Xi X_(i+W) + Yi Y_(i+1) +Yi Y_(i+W)]
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ [S+_i S+_(i+1) + S-_iS-_(i+1) + S+_i S+_(i+W) + S-_iS-_(i+W)]
"""
function XYmodel_2D_helical(wid::Integer = 1; expand::Bool = false)
    phy_dim = 2; vir_dim = 2*wid
    sp = pauli(PPlus); sm = pauli(PMinus);
    Tp = Ising_CMPO(0.5, sp, sp, wid)
    Tm = Ising_CMPO(0.5, sm, sm, wid)
    Tmatrix = cat(Tp, Tm)
    Ut = nothing
    if expand
        Tmatrix = Tmatrix |> expand_cmpo
        Ut, vir_dim = generalUt(vir_dim)
    end
    return PhysModel(Tmatrix, phy_dim, vir_dim, Ut)
end


"""
2D Heisenberg XXZ model helical boundary condition
    H = ∑ [Xi X_(i+1) +Xi X_(i+W) + Yi Y_(i+1) +Yi Y_(i+W)] + Δ [Zi Z_(i+1) +Zi Z_(i+W)]
after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5 ∑ [S+_i S+_(i+1) + S-_iS-_(i+1) + S+_i S+_(i+W) + S-_iS-_(i+W)] - Δ [Zi Z_(i+1) +Zi Z_(i+W)]
"""
function XXZmodel_2D_helical(Δ::Real, wid::Integer = 1; expand::Bool=false)
    phy_dim = 2; vir_dim = 3*wid
    sp = pauli(PPlus); sm = pauli(PMinus); sz = 0.5 * pauli(PZ)
    Tp = Ising_CMPO(0.5, sp, sp, wid)
    Tm = Ising_CMPO(0.5, sm, sm, wid)
    Tz = Ising_CMPO(Δ, sz, sz, wid)
    Tmatrix = cat(cat(Tp, Tm), Tz) 
    Ut = nothing
    if expand 
        Tmatrix = Tmatrix |> expand_cmpo
        Ut, vir_dim = generalUt(vir_dim)
    end
    return PhysModel(Tmatrix, phy_dim, vir_dim, Ut)
end


"""
2D general XXZ model helical boundary condition
    H = ∑ Jxy[Xi X_(i+1) +Xi X_(i+W) + Yi Y_(i+1) +Yi Y_(i+W)] + Jz [Zi Z_(i+1) +Zi Z_(i+W)]
after unitary transformation : U=exp(iπSy) on odd sites:
    H = 0.5Jxy ∑ [S+_i S-_(i+1) + S-_iS+_(i+1) + S+_i S-_(i+W) + S-_iS+_(i+W)] + Δ [Zi Z_(i+1) +Zi Z_(i+W)]
"""
function XXZmodel_2D_helical(Jxy::Real, Jz::Real, wid::Integer = 1; expand::Bool=false)
    phy_dim = 2; vir_dim = 3*wid
    sp = pauli(PPlus); sm = pauli(PMinus); sz = 0.5 * pauli(PZ)
    Tp = Ising_CMPO(-0.5*Jxy, sp, sm, wid)
    Tm = Ising_CMPO(-0.5*Jxy, sm, sp, wid)
    Tz = Ising_CMPO(-Jz, sz, sz, wid)
    Tmatrix = cat(cat(Tp, Tm), Tz) 
    Ut = nothing
    if expand 
        Tmatrix = Tmatrix |> expand_cmpo
        Ut, vir_dim = generalUt(vir_dim)
    end
    return PhysModel(Tmatrix, phy_dim, vir_dim, Ut)
end
#end module PhysicalModels
