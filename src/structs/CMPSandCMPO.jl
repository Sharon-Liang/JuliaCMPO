#module structs
"""
    AbstractCTensor 

Abstract type of CMPS and CMPO tensors
"""
abstract type AbstractCTensor end


"""
    AbstractCMPS{T, S, U} <: AbstractCTensor

Abstract type of all CMPS, the structure of a CMPS tensor:
--       --
| I + ϵQ  |
|         |
|   |     |
| √ϵR     |
|   |     |
--       --
"""
abstract type AbstractCMPS{T, S, U} <: AbstractCTensor end


"""
    AbstractCMPO{T,S,U,V} <: AbstractCTensor

Abstract type of all CMPO, the structure of a CMPO tensor:
--                 --
| I + ϵQ  -- √ϵL -- |
|                   |
|   |               |
| √ϵR        P      |
|   |               |
--                 --
"""
abstract type AbstractCMPO{T,S,U,V} <: AbstractCTensor end


""" 
    CMPS <: AbstractCMPS
"""
@with_kw struct CMPS{T<:Number, S<:Matrix, U<:Array} <: AbstractCMPS{T, S, U}
    Q::S # Dimension: χ × χ 
    R::U # Dimension: χ × χ × vir_dim
    CMPS{T,S,U}(Q::Matrix{T}, R::Array{T}) where {T,S,U} = new(Q,R)
end
CMPS(Q::Matrix{T}, R::Array{T}) where {T} = CMPS{T,typeof(Q), typeof(R)}(Q,R)


"""
    CuCMPS <: AbstractCMPS
"""
@with_kw struct CuCMPS{T<:Number, S<:CuMatrix, U<:CuArray} <: AbstractCMPS{T, S, U}
    Q::S # Dimension: χ × χ 
    R::U# Dimension: χ × χ × vir_dim
    CuCMPS{T,S,U}(Q::CuMatrix{T}, R::CuArray{T}) where {T,S,U} = new(Q,R)
end
CuCMPS(Q::CuMatrix{T}, R::CuArray{T}) where {T} = CuCMPS{T, typeof(Q), typeof(R)}(Q,R)


"""
    CMPO <: AbstractCMPO
"""
@with_kw struct CMPO{T<:Number, S<:Matrix, U<:Array, V<:Array} <: AbstractCMPO{T,S,U,V}
    Q::S # Dimension: phy_dim × phy_dim
    R::U # Column : phy_dim × phy_dim × vir_dim
    L::U  # Row: phy_dim × phy_dim × vir_dim
    P::V  # long-range interaction: phy_dim × phy_dim × vir_dim × vir_dim
    CMPO{T,S,U,V}(Q::Matrix{T}, R::Array{T}, L::Array{T}, P::Array{T}) where {T,S,U,V} = new(Q,R,L,P)
end
CMPO(Q::Matrix{T}, R::Array{T}, L::Array{T}, P::Array{T}) where {T} = 
    CMPO{T, typeof(Q), typeof(R), typeof(P)}(Q,R,L,P)

"""
    CuCMPO <: AbstractCMPO
"""
@with_kw struct CuCMPO{T<:Number, S<:CuMatrix, U<:CuArray, V<:CuArray} <: AbstractCMPO{T,S,U,V}
    Q::S # Dimension: phy_dim × phy_dim
    R::U # Column : phy_dim × phy_dim × vir_dim
    L::U  # Row: phy_dim × phy_dim × vir_dim
    P::V  # long-range interaction: phy_dim × phy_dim × vir_dim × vir_dim
    CuCMPO{T,S,U,V}(Q::CuMatrix{T}, R::CuArray{T}, L::CuArray{T}, P::CuArray{T}) where {T,S,U,V} = 
        new(Q,R,L,P)
end
CuCMPO(Q::CuMatrix{T}, R::CuArray{T}, L::CuArray{T}, P::CuArray{T}) where {T} = 
        CuCMPO{T, typeof(Q), typeof(R), typeof(P)}(Q,R,L,P)


"""
    CMPS/CMPO generation function
"""
cmps_generate(QR::Array...) = CMPS(QR...)
cmps_generate(QR::CuArray...) = CuCMPS(QR...)
cmps_generate(;Q, R) = cmps_generate(Q, R)

cmpo_generate(QRLP::Array...) = CMPO(QRLP...)
cmpo_generate(QRLP::CuArray...) = CuCMPO(QRLP...)
cmpo_generate(;Q, R, L, P) = cmps_generate(Q, R, L, P)


"""
    CTensor(x::AbstractCTensor)

Convert AbstractCTensors to `CMPS`/`CMPO`
"""
CTensor(x::Union{CMPS, CMPO}) = x
CTensor(x::CuCMPS) = CMPS(Matrix(x.Q), Array(x.R))
CTensor(x::CuCMPO) = CMPO(Matrix(x.Q), Array(x.R), Array(x.L), Array(x.P))

"""
    CuCTensor(x::AbstractCTensor)

Convert AbstractCTensors to `CuCMPS`/`CuCMPO`
"""
CuCTensor(x::Union{CuCMPS, CuCMPO}) = x
CuCTensor(x::CMPS) = CuCMPS(CuMatrix(x.Q), CuArray(x.R))
CuCTensor(x::CMPO) = CuCMPO(CuMatrix(x.Q), CuArray(x.R), CuArray(x.L), CuArray(x.P))


bond_dimension(x::AbstractCTensor) = size(x.Q, 1)
virtual_bond_dimension(x::AbstractCTensor) = Integer(length(x.R)/length(x.Q)) + 1

#end #module structs

