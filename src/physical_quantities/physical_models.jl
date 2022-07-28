#module PhysicalModels
"""
    Ising_cmpo(J, ol, or[, wid = 1])

Generate ising type CMPO ``H = -J ôl ôr`` block.
"""
function ising_cmpo(J::Real, ol::AbstractArray, or::AbstractArray, wid::Integer = 1)
    sgn = sign(J); val =√(abs(J))
    o2 = zeros(eltype(ol), 2, 2)
    i2 = oneunit(o2)
    #construct L and R
    ol = val .* ol; or = sgn .* val .* or
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
        
    return CMPO(o2, R, L, P)
end


"""
    cat(o1, o2)

Cat two CMPO blocks.
"""
function Base.cat(o1::AbstractCMPO, o2::AbstractCMPO)
    pd = size(o1.Q,1)
    typeof(o1.P) <:AbstractMatrix ? vd1 = 1 : vd1 = size(o1.P,3)
    typeof(o2.P) <:AbstractMatrix ? vd2 = 1 : vd2 = size(o2.P,3)

    Q = zeros(eltype(o1.Q), pd, pd)
    R = cat(o1.R, o2.R, dims = 3)
    L = cat(o1.L, o2.L, dims = 3)

    pl = cat(o1.P, zeros(eltype(o1.P), pd, pd, vd2, vd1), dims = 3)
    pr = cat(zeros(eltype(o2.P), pd, pd, vd1, vd2), o2.P, dims = 3)
    P = cat(pl, pr, dims = 4)

    return CMPO(Q, R, L, P)
end


"""
    TFIsing(J, Γ[; testfield = nothing, η=1.e-2])

Generate the CMPO of 1D nearest neighbor transvers field Ising model:
```math
H = -∑ J Zi Zj - ∑ Γ Xi
```
where ``X,Y,Z`` are Pauli matrices.
"""
function TFIsing(J::Real, Γ::Real; testfield = nothing, η::Float64 = 1.e-2)
    testfield === nothing ? h = zeros(2,2) : h = η .* pauli(testfield)
    cmpo = CMPO(Γ*pauli(PX)+h, √J*pauli(PZ), √J*pauli(PZ), zeros(2,2))
    return cmpo
end


"""
    XYmodel()

Generate the CMPO of 1D nearest neighbor XY model
```math
H = ∑ Xi Xj + Yi Yj
```
after unitary transformation ``U=exp(iπSy)`` on odd sites:
```math
H = -0.5 ∑ (S₊ S₊ + S₋ S₋)
```
"""
function XYmodel()
    sp = pauli(PPlus); sm = pauli(PMinus);
    L = zeros(2, 2, 2)
    L[:,:,1] = 1/√2 * sp ; L[:,:,2] = 1/√2 * sm
    R = zeros(2, 2, 2)
    R[:,:,1] = 1/√2 * sp ; R[:,:,2] = 1/√2 * sm
    Q = zeros(2, 2)
    P = zeros(2, 2, 2, 2)
    return CMPO(Q, R, L, P)
end


"""
    XXZmodel(Δ)

Generate the CMPO of 1D Heisenberg XXZ model
```math
H = ∑ Xi Xj + Yi Yj + Δ Zi Zj
```
after unitary transformation ``U=exp(iπSy)`` on odd sites:
```math
H = -0.5 ∑ (S₊ S₊ + S₋ S₋) - ∑ Δ Zi Zj
```
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
        return CMPO(Q,R,L,P)
    end
end


"""
    TFIsing_2D_helical(J, Γ[, wid = 1; testfield = nothing, η=1.e-2])

Generate the CMPO of nearest neighbor transvers field Ising model
on a 2D cylinder of width `W` with helical boundary condition:
```math
H = -∑ J [Zi Z_(i+1) +Zi Z_(i+W)]  - ∑ Γ Xi
```
where ``X,Y,Z`` are Pauli matrices.
"""
function TFIsing_2D_helical(J::Real, Γ::Real, wid::Integer = 1; testfield = nothing, η::Float64 = 1.e-2)
    testfield === nothing ? h = zeros(2,2) : h = η .* pauli(field)
    Q = Γ*pauli(PX)+h
    Tm = ising_cmpo(J, pauli(PZ), pauli(PZ), wid)
    return CMPO(Q, Tm.R, Tm.L, Tm.P)
end


"""
    XYmodel_2D_helical([wid = 1])

Generate the CMPO of nearest neighbor XY model
on a 2D cylinder of width `W` with helical boundary condition:
```math
H = -0.5 ∑ [S₊(i) S₊(i+1) + S₋(i)S₋(i+1) + S₊(i)S₊(i+W) + S₋(i)S₋(i+W)]
```    
"""
function XYmodel_2D_helical(wid::Integer = 1)
    sp = pauli(PPlus); sm = pauli(PMinus);
    Tp = ising_cmpo(0.5, sp, sp, wid)
    Tm = ising_cmpo(0.5, sm, sm, wid)
    return cat(Tp, Tm)
end


"""
    XXZmodel_2D_helical(Δ[, wid = 1])

Generate the CMPO of nearest neighbor Heisenberg XXZ model
on a 2D cylinder of width `W` with helical boundary condition:
```math 
H = -0.5 ∑ [S₊(i) S₊(i+1) + S₋(i)S₋(i+1) + S₊(i)S₊(i+W) + S₋(i)S₋(i+W)] - Δ [Z(i) Z(i+1) + Z(i) Z(i+W)]
```
"""
function XXZmodel_2D_helical(Δ::Real, wid::Integer = 1)
    sp = pauli(PPlus); sm = pauli(PMinus); sz = 0.5 * pauli(PZ)
    Tp = ising_cmpo(0.5, sp, sp, wid)
    Tm = ising_cmpo(0.5, sm, sm, wid)
    Tz = ising_cmpo(Δ, sz, sz, wid)
    return cat(cat(Tp, Tm), Tz) 
end


"""
    XXZmodel_2D_helical(Jxy, Jz[, wid = 1])

Generate the CMPO of nearest neighbor Heisenberg XXZ model
on a 2D cylinder of width `W` with helical boundary condition:
```math
H = -0.5 Jxy ∑ [S₊(i) S₊(i+1) + S₋(i)S₋(i+1) + S₊(i)S₊(i+W) + S₋(i)S₋(i+W)] - Jz [Z(i) Z(i+1) + Z(i) Z(i+W)]
```
"""
function XXZmodel_2D_helical(Jxy::Real, Jz::Real, wid::Integer = 1)
    sp = pauli(PPlus); sm = pauli(PMinus); sz = 0.5 * pauli(PZ)
    Tp = ising_cmpo(-0.5*Jxy, sp, sm, wid)
    Tm = ising_cmpo(-0.5*Jxy, sm, sp, wid)
    Tz = ising_cmpo(-Jz, sz, sz, wid)
    return cat(cat(Tp, Tm), Tz) 
end
#end module PhysicalModels
