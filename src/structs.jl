#module structs
"""
χ::Int64       # χ: bond dimension along imaginary time direction
vir_dim::Int64 # vir_dim+1: virtual bond dimension of the horizontal legs
phy_dim::Int64 # phy_dim: bond dimension of physical legs
"""
struct CMPS{T<:Number}
    Q::Matrix{T} # Dimension: χ × χ 
    R::Array{T}  # Dimension: χ × χ × vir_dim
end

struct CMPO{T<:Number}
    Q::Matrix{T} # Dimension: phy_dim × phy_dim
    R::Array{T}  # Column : phy_dim × phy_dim × vir_dim
    L::Array{T}  # Row: phy_dim × phy_dim × vir_dim
    P::Array{T}  # long-range interaction: phy_dim × phy_dim × vir_dim × vir_dim
end

struct PhysModel
    Tmatrix::CMPO # local Transfer matrix
    phy_dim::Int64  # phy_dim: bond dimension of physical legs
    vir_dim::Int64  # vir_dim: virtual bond dimension of the horizontal legs
    Ut::Matrix{T} where T <: Number # unitrary transformation: Ut^+ Tm Ut = transpose(Tm)
end

struct MeraUpdateStep{Ti<:Integer, T<:Real, Tf<:Real}
    SN::Ti
    θ::T
    loss_diff::Tf
    fidelity::Tf
end

MeraUpdateTrace = Vector{MeraUpdateStep}

#end module structs
