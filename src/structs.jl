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
abstract type AbstractCMPS{T} <: AbstractCTensor end
abstract type AbstractCMPO{T} <: AbstractCTensor end
struct CMPS{T<:Number} <: AbstractCMPS{T}
    Q::Matrix{T} # Dimension: χ × χ 
    R::Array{T} # Dimension: χ × χ × vir_dim
end

struct CuCMPS{T<:Number} <: AbstractCMPS{T}
    Q::CuMatrix{T} # Dimension: χ × χ 
    R::CuArray{T} # Dimension: χ × χ × vir_dim
end


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
struct CMPO{T<:Number} <: AbstractCMPO{T}
    Q::Matrix{T} # Dimension: phy_dim × phy_dim
    R::Array{T} # Column : phy_dim × phy_dim × vir_dim
    L::Array{T}  # Row: phy_dim × phy_dim × vir_dim
    P::Array{T}  # long-range interaction: phy_dim × phy_dim × vir_dim × vir_dim
end

struct CuCMPO{T<:Number} <: AbstractCMPO{T}
    Q::CuMatrix{T} # Dimension: phy_dim × phy_dim
    R::CuArray{T} # Column : phy_dim × phy_dim × vir_dim
    L::CuArray{T}  # Row: phy_dim × phy_dim × vir_dim
    P::CuArray{T}  # long-range interaction: phy_dim × phy_dim × vir_dim × vir_dim
end


"""
    CMPS/CMPO generation function
"""
CMPS_generate(QR::Array{T}...) where T<:Number = CMPS(QR...)
CMPS_generate(QR::CuArray{T}...) where T<:Number = CuCMPS(QR...)

CMPO_generate(QRLP::Array{T}...) where T<:Number = CMPO(QRLP...)
CMPO_generate(QRLP::CuArray{T}...) where T<:Number = CuCMPO(QRLP...)

"""
    `CTensor(x)`: convert CTensors to CMPS/CMPO
    `CuCTensor(x)` : convert CTensors to CuCMPS/CuCMPO
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

#end #module structs
