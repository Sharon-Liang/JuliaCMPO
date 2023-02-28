"""
    AbstractModel

An abstract type representing the physical models
"""
abstract type AbstractModel end


"""
    model(a::AbstractModel)

Constract model cMPO.
"""
model(a::AbstractModel) = "Undefined."


#=
### *1D Models*
=#
"""
    TFIsingChain

1D transverse field Ising model with nearest neighbor interactions:
```math
    H = - J ∑ σᶻᵢ σᶻⱼ - Γ ∑ σˣᵢ
```
"""
@with_kw struct TFIsingChain{T<:Real} <: AbstractModel 
    J::T 
    Γ::T 
end

function model(M::TFIsingChain)
    @unpack J, Γ = M
    Q = Γ*pauli(PX)
    R = √J*pauli(PZ)
    L = √J*pauli(PZ)
    P = zeros(2,2)
    return CMPO(Q, R, L, P)
end


"""
    XXChain

1D nearst neighbor XX model:
```math
    H = ∑ (SˣᵢSˣⱼ + Sʸᵢ Sʸⱼ)
```
after unitary transformation ``U = exp(iπSʸ)`` on odd sites, the Hamiltoinan becomes:
```math
    H = -0.5 ∑ (S⁺ᵢS⁺ⱼ + S⁻ᵢ S⁻ⱼ)
```
We implement the later form in this package.
"""
struct XXChain <: AbstractModel end

function model(::XXChain)
    sp = pauli(PPlus) 
    sm = pauli(PMinus)
    Q = zeros(2,2)
    R = cat(1/√2 * sp, 1/√2 * sm, dims=3)
    L = cat(1/√2 * sp, 1/√2 * sm, dims=3)
    P = zeros(2, 2, 2, 2)
    return CMPO(Q, R, L, P)
end


"""
    XXZChain

1D nearst neighbor Heisenberg XXZ model:
```math
    H = ∑ (SˣᵢSˣⱼ + Sʸᵢ Sʸⱼ + Δ Sᶻᵢ Sᶻⱼ)
```
after unitary transformation ``U = exp(iπSʸ)`` on odd sites, the Hamiltoinan becomes:
```math
    H = -0.5 ∑ (S⁺ᵢS⁺ⱼ + S⁻ᵢ S⁻ⱼ) - Δ ∑ Sᶻᵢ Sᶻⱼ
```
We implement the later form in this package.
"""
@with_kw struct XXZChain{T<:Real} <: AbstractModel 
    Δ::T
end

function model(M::XXZChain) 
    Δ = M.Δ  
    if Δ == 0
        return model(XXmode())
    else
        sgn = sign(Δ); val = √abs(Δ)
        sp = pauli(PPlus)
        sm = pauli(PMinus)
        sz = 0.5 * pauli(PZ) 
        Q = zeros(2, 2)
        R = cat(1/√2 * sp, 1/√2 * sm, val * sz, dims=3)
        L = cat(1/√2 * sp, 1/√2 * sm, sgn * val * sz, dims=3)
        P = zeros(2, 2, 3, 3)
        return CMPO(Q, R, L, P)
    end
end


"""
    J1J2Chain

1D `J₁-J₂` model, where `J₁` is the nearest neighbor interaction while `J₂` is the next nearest neighbor interaction, the Hamiltoinan is:
```math
    H = J₁ ∑ (SˣᵢSˣⱼ - iSʸᵢ iSʸⱼ + Sᶻᵢ Sᶻⱼ)
            + J₂ ∑ (SˣᵢSˣₖ - iSʸᵢ iSʸₖ + Sᶻᵢ Sᶻₖ)
```
where `j = i+1` and `k=i+2`.
"""
@with_kw struct J1J2Chain{T<:Real} <: AbstractModel 
    J1::T
    J2::T
end


"""
    _j1_j2_block(J1::Real, J2::Real, p::PauliMatrixName)

Construct CMPO blocks for terms related to operator `pauli(p)`, e.g.
```math
    H = J₁ ∑ SˣᵢSˣⱼ + J₂ ∑ SˣᵢSˣₖ 
```
where `j = i+1` and `k=i+2`.
"""
function _j1_j2_block(J1, J2, p::PauliMatrixName)
    op = 0.5 * pauli(p)
    s1 = sign(J1); v1 =√(abs(J1))
    s2 = sign(J2); v2 =√(abs(J2))

    Q = zeros(2, 2)
    R = cat( v1 * op, zeros(2,2), dims=3)
    L = cat(-s1 * v1 * op, -s2 * v2 * op, dims=3)
    P = zeros(2,2,2,2); P[:,:,2,1] = v2/v1 * oneunit(op)
    return CMPO(Q, R, L, P)
end


function model(M::J1J2Chain)
    @unpack J1, J2 = M
    Tx = _j1_j2_block(J1, J2, PX)
    Ty = _j1_j2_block(-J1, -J2, iPY)
    Tz = _j1_j2_block(J1, J2, PZ)
    return cat(Tx, Ty, Tz)
end



#=
### *2D Models* : *Square Lattice, Helical boundary condition*
=#
"""
    _ising_2D_block(J::Real, ol::AbstractArray, or::AbstractArray[, W::Int = 1])

Construct CMPO blocks of Hamiltonians with nearst neighbor Ising type interactions on a 2D square lattice cylinder, the Wth of the cylinder is `W`, and `W=1` corresponding to a Ising Chain. The corresponding quansi 1D Hamiltoinan is:
```math
    H = - J ∑ σᶻᵢ σᶻⱼ - J ∑ σᶻᵢ σᶻₖ
```
where `j = i+1` and `k=i+W`.
"""
function _ising_2D_block(J::Real, ol::AbstractArray, or::AbstractArray, W::Int = 1)
    sgn = sign(J); val =√(abs(J))
    o2 = zeros(eltype(ol), size(ol))
    i2 = oneunit(o2)

    #construct L and R
    ol = val .* ol
    or = sgn .* val .* or

    L = ol
    R = or
    for i = 2: W
        R = cat(R, o2, dims = 3)
        i == W ? L = cat(L, ol, dims=3) : L = cat(L, o2, dims = 3)
    end

    #construct P 
    P = []
    for i = 1: W
        pcol = o2
        for j = 2: W
            j == i+1 ? pcol = cat(pcol, i2, dims = 3) : pcol = cat(pcol, o2, dims = 3)
        end
        i == 1 ? P = pcol : P = cat(P, pcol, dims = 4)
    end
        
    return CMPO(o2, R, L, P)
end


"""
    TFIsingSquareHelical

2D nearst neighbor transverve field Ising model on a square lattice that wrapped on a cylinder. The quansi 1D Hamiltoinan is:
```math
    H = -J ∑ ( σᶻᵢ σᶻⱼ - σᶻᵢ σᶻₖ) - Γ ∑ σˣᵢ
```
where `j = i+1` and `k=i+W`.
"""
struct TFIsingSquareHelical{T<:Real} <: AbstractModel 
    J::T
    Γ::T 
    W::Int
end

function model(M::TFIsingSquareHelical)
    @unpack J, Γ, W = M
    Q = Γ * pauli(PX)
    T = _ising_2D_block(J, pauli(PZ), pauli(PZ), W)
    return CMPO(Q, T.R, T.L, T.P)
end 


"""
    XXSquareHelical

2D nearst neighbor XX model on a square lattice that wrapped on a cylinder. The quansi 1D Hamiltoinan is:
```math
    H = -0.5 ∑ (S⁺ᵢS⁺ⱼ + S⁻ᵢ S⁻ⱼ) -0.5 ∑ (S⁺ᵢS⁺ₖ + S⁻ᵢ S⁻ₖ)
```
where `j = i+1` and `k=i+W`.
"""
struct XXSquareHelical <: AbstractModel 
    W::Int
end

function model(M::XXSquareHelical)
    W = M.W
    sp = pauli(PPlus)
    sm = pauli(PMinus)
    Tp = _ising_2D_block(0.5, sp, sp, W)
    Tm = _ising_2D_block(0.5, sm, sm, W)
    return cat(Tp, Tm)
end


"""
    XXZSquareHelical

2D nearst neighbor XXZ model on a square lattice that wrapped on a cylinder. The quansi 1D Hamiltoinan is:
```math
    H = -0.5 ∑ (S⁺ᵢS⁺ⱼ + S⁻ᵢ S⁻ⱼ) - Δ ∑ Sᶻᵢ Sᶻⱼ
        -0.5 ∑ (S⁺ᵢS⁺ₖ + S⁻ᵢ S⁻ₖ) - Δ ∑ Sᶻᵢ Sᶻₖ
```
where `j = i+1` and `k=i+W`.
"""
struct XXZSquareHelical{T<:Real} <: AbstractModel 
    Δ::T
    W::Int
end

function model(M::XXZSquareHelical)
    @unpack Δ, W = M
    sp = pauli(PPlus)
    sm = pauli(PMinus)
    sz = 0.5 * pauli(PZ)

    Tp = _ising_2D_block(0.5, sp, sp, W)
    Tm = _ising_2D_block(0.5, sm, sm, W)
    Tz = _ising_2D_block(Δ, sz, sz, W)
    return cat(Tp, Tm, Tz)
end