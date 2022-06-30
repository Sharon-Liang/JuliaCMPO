#module structs
""" `CMPS`: the struct for cMPS
    `vir_dim::Int64` : vir_dim+1 is the virtual bond dimension of the horizontal legs
    `χ::Int64` : bond dimension along imaginary time direction
the structure of cMPO 
    --       --
    | I + ϵQ  |
    |         |
    |   |     |
    | √ϵR     |
    |   |     |
    --       --
"""
abstract type AbstractCTensor end
abstract type AbstractCMPS{T, S, U} <: AbstractCTensor end
abstract type AbstractCMPO{T, S, U, V} <: AbstractCTensor end
@with_kw struct CMPS{T<:Number, S<:Matrix, U<:Array} <: AbstractCMPS{T, S, U}
    Q::S # Dimension: χ × χ 
    R::U # Dimension: χ × χ × vir_dim
    CMPS{T,S,U}(Q::Matrix{T}, R::Array{T}) where {T,S,U} = new(Q,R)
end
CMPS(Q::Matrix{T}, R::Array{T}) where {T} = CMPS{T,typeof(Q), typeof(R)}(Q,R)

@with_kw struct CuCMPS{T<:Number, S<:CuMatrix, U<:CuArray} <: AbstractCMPS{T, S, U}
    Q::S # Dimension: χ × χ 
    R::U# Dimension: χ × χ × vir_dim
    CuCMPS{T,S,U}(Q::CuMatrix{T}, R::CuArray{T}) where {T,S,U} = new(Q,R)
end
CuCMPS(Q::CuMatrix{T}, R::CuArray{T}) where {T} = CuCMPS{T, typeof(Q), typeof(R)}(Q,R)


""" `CMPO`: the struct for cMPO
    `vir_dim::Int64` : vir_dim+1 is the virtual bond dimension of the horizontal legs
    `phy_dim::Int64` : bond dimension of physical legs
the structure of cMPO 
    --                 --
    | I + ϵQ  -- √ϵL -- |
    |                   |
    |   |               |
    | √ϵR        P      |
    |   |               |
    --                 --
"""
@with_kw struct CMPO{T<:Number, S<:Matrix, U<:Array, V<:Array} <: AbstractCMPO{T, S, U, V}
    Q::S # Dimension: phy_dim × phy_dim
    R::U # Column : phy_dim × phy_dim × vir_dim
    L::U  # Row: phy_dim × phy_dim × vir_dim
    P::V  # long-range interaction: phy_dim × phy_dim × vir_dim × vir_dim
    CMPO{T,S,U,V}(Q::Matrix{T}, R::Array{T}, L::Array{T}, P::Array{T}) where {T,S,U,V} = new(Q,R,L,P)
end
CMPO(Q::Matrix{T}, R::Array{T}, L::Array{T}, P::Array{T}) where {T} = 
    CMPO{T, typeof(Q), typeof(R), typeof(P)}(Q,R,L,P)

@with_kw struct CuCMPO{T<:Number, S<:CuMatrix, U<:CuArray, V<:CuArray} <: AbstractCMPO{T, S, U, V}
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
CMPS_generate(QR::Array{T}...) where T<:Number = CMPS(QR...)
CMPS_generate(QR::CuArray{T}...) where T<:Number = CuCMPS(QR...)
CMPS_generate(;Q, R) = CMPS_generate(Q, R)

CMPO_generate(QRLP::Array{T}...) where T<:Number = CMPO(QRLP...)
CMPO_generate(QRLP::CuArray{T}...) where T<:Number = CuCMPO(QRLP...)
CMPO_generate(;Q, R, L, P) = CMPS_generate(Q, R, L, P)


"""
    `CTensor(x)`: convert CTensors to CMPS/CMPO
    `CuCTensor(x)` : convert CTensors to CuCMPS/CuCMPOs
"""
CTensor(x::Union{CMPS, CMPO}) = x
function CTensor(x::T) where T<:AbstractCTensor
    fields = fieldnames(T)
    args =  Tuple(convert(Array, getfield(x, field)) for field in fields)
    length(fields) == 2 ? CMPS(args...) : CMPO(args...)
end

CuCTensor(x::Union{CuCMPS, CuCMPO}) = x
function CuCTensor(x::T) where T<:AbstractCTensor
    fields = fieldnames(T)
    args =  Tuple(convert(CuArray, getfield(x, field)) for field in fields)
    length(fields) == 2 ? CuCMPS(args...) : CuCMPO(args...)
end

bond_dimension(x::AbstractCTensor) = size(x.Q)[1]
virtual_dimension(x::AbstractCTensor) = Integer(length(x.R)/length(x.Q)) + 1
#end #module structs

