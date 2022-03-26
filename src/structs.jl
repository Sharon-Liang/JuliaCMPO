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

struct PhysModel{T}
    Tmatrix::CMPO # local Transfer matrix
    phy_dim::Int64  # phy_dim: bond dimension of physical legs
    vir_dim::Int64  # vir_dim: virtual bond dimension of the horizontal legs
    Ut::T # unitrary transformation: Ut^+ Tm Ut = transpose(Tm)
end

struct MeraUpdateStep{Ti<:Integer, T<:Real, Tf<:Real}
    SN::Ti
    θ::T
    loss_diff::Tf
    fidelity::Tf
end

MeraUpdateTrace = Vector{MeraUpdateStep}

struct MeraUpdateResult
    ψ::CMPS
    trace::MeraUpdateTrace
end


"""
atol: tolerance of absolute difference of fidelity(ψ, ψ0, β, Normalize = true) and 1.0
ldiff_tol: tolerance of absolute difference of the value of loss function between two MERA update steps
"""
struct MeraUpdateOptions{T<:Number}
    atol::T
    ldiff_tol::T
    maxiter::Int64
    interpolate::Bool
    store_trace::Bool
    show_trace::Bool
end

function MeraUpdateOptions(;
    atol = 1.e-5,
    ldiff_tol = 1.e-12,
    maxiter = 100,
    interpolate = true,
    store_trace = false,
    show_trace = false)
    MeraUpdateOptions(atol, ldiff_tol, maxiter, interpolate,
        store_trace, show_trace)
end

struct CompressResult{T<:Number, Tf}
    ψ::CMPS
    fidelity_initial::T
    fidelity_final::T
    optim_result::Tf
end

#end module structs
