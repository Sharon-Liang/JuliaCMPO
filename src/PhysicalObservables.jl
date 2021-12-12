#module PhysicalObservables
"""
    make operator ⊢o⊣
"""
function make_operator(Op::AbstractMatrix, dim::Int64)
    T = eltype(Op)
    eye = Matrix{T}(I, dim, dim)
    return eye ⊗ Op ⊗ eye
end

function make_operator(Op::AbstractMatrix, ψ::CMPS)
    T = eltype(Op)
    eye = Matrix{T}(I, size(ψ.Q))
    return eye ⊗ Op ⊗ eye
end


"""
    The thermal average of local opeartors ⊢o⊣ with respect to K = ψ * W * ψ
"""
function thermal_average(Op::AbstractMatrix, ψ::CMPS, W::CMPO, β::Real)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(-β*K)
    m = maximum(e); e = e .- m
    Op = v' * Op * v
    den = exp.(e) |> sum
    num = exp.(e) .* diag(Op) |> sum
    return num/den
end

"""
    The thermal average of local opeartors ⊢o⊣ with respect to K = ψ * ψ
"""
function thermal_average(Op::AbstractMatrix, ψ::CMPS, β::Real)
    K = ψ * ψ |> symmetrize |> Hermitian
    e, v = eigen(-β*K)
    m = maximum(e) ; e = e .- m
    Op = v' * Op * v
    den = exp.(e) |> sum
    num = exp.(e) .* diag(Op) |> sum
    return num/den
end


"""
    Partitian function
"""
function partitian(ψ::CMPS, W::CMPO, β::Real)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    H = ψ * ψ |> symmetrize |> Hermitian
    num = trexp(-β*K)
    den = trexp(-β*H)
    return exp(num.max - den.max) * num.res/den.res
end


"""
    Free energy
"""
function free_energy(ψ::CMPS, W::CMPO, β::Real)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    H = ψ * ψ |> symmetrize |> Hermitian
    res = logtrexp(-β*K)- logtrexp(-β*H)
    return -1/β * res
end

function free_energy(param::Array{<:Number,3}, W::CMPO, β::Real)
    free_energy(tocmps(param), W, β)
end

function free_energy(param::Vector{<:Number}, dim::Tuple, W::CMPO, β::Real)
    free_energy(tocmps(param, dim), W, β)
end

function energy(ψ::CMPS, W::CMPO, β::Real)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    H = ψ * ψ |> symmetrize |> Hermitian
    res = thermal_average(K, ψ, W, β) - thermal_average(H, ψ, β)
    return res
end

function specific_heat(ψ::CMPS, W::CMPO, β::Real; method::Symbol = :ndiff)
    if method == :adiff
        K = ψ * W * ψ |> symmetrize |> Hermitian
        H = ψ * ψ |> symmetrize |> Hermitian
        K2 = K * K
        H2 = H * H
        c = thermal_average(K2, ψ, W, β) - thermal_average(K, ψ, W, β)^2
        c -= thermal_average(H2, ψ, β) - thermal_average(H, ψ, β)^2
    elseif method == :ndiff
        e = b -> energy(ψ, W, b)
        c = -central_fdm(5, 1)(e, β)
    else @error "method should be :adiff or :ndiff"
    end
    return β^2 * c
end


function entropy(ψ::CMPS, W::CMPO, β::Real)
    s = energy(ψ,W,β) - free_energy(ψ,W,β)
    return β*s
end


"""
The local two-time correlation functions
"""
function correlation_2time(τ::Number, A::AbstractMatrix,B::AbstractMatrix,
                           ψ::CMPS, W::CMPO, β::Real)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    #m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    #den = exp.(-β * e .- m) |> sum
    den = exp.(-β * e) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        #num += exp(-β*e[i]- m + τ*(e[i] - e[j])) * A[i,j] * B[j,i]
        num += exp(-β*e[i] + τ*(e[i] - e[j])) * A[i,j] * B[j,i]
    end
    return num/den
end


"""
check anomalous term of bosonic Masubara correlations
n = 0  and Em = En terms
"""
function check_anomalous_term(A::AbstractMatrix,B::AbstractMatrix,
    ψ::CMPS, W::CMPO, β::Real)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    #m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    #den = exp.(-β * e .- m) |> sum
    den = exp.(-β * e) |> sum
    c = 0.
    for i = 1: length(e), j = 1: length(e)
        if e[i] == e[j]
            c += exp(-β*e[j]) * A[i,j] * B[j,i]
        end
    end
    return β*c/den       
end


"""
A useful function: f = e^(-b*e1) - e^(-b*e2) / (e2 - e1)
"""
function diffaddexp(b::Real, e1::Real, e2::Real)
    if abs(e2 - e1) < 1.e-10
        return exp(-b*e1) * b
    else
        num = exp(-b*e1) - exp(-b*e2)
        den = e2 - e1
        return num/den
    end
end

"""
Masubara frequency Green's functions: defalt type = :b
"""
function Masubara_freq_GF(n::Integer, A::AbstractMatrix,B::AbstractMatrix,
                        ψ::CMPS, W::CMPO, β::Real)
    λ = 1.0
    ωn = Masubara_freq(n,β,type=:b)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    #m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    #den = exp.(-β * e .- m) |> sum
    num = 0.0
    if ωn != 0
        for i = 1: length(e), j = 1: length(e)
            up = exp(-β*e[i]) - λ*exp(-β*e[j])
            #up = exp(-β*e[i]-m) - λ*exp(-β*e[j]-m)
            up = up * A[i,j] * B[j,i]
            down = 1.0im * ωn - e[j] + e[i]
            num += up/down
        end
    else
        for i = 1: length(e), j = 1: length(e)
            num -= A[i,j]*B[j,i]*diffaddexp(β,e[i],e[j])
        end
    end
    return num/den
end

function Masubara_freq_GF(z::ComplexF64, A::AbstractMatrix,B::AbstractMatrix,
    ψ::CMPS, W::CMPO, β::Real)
    λ = 1.0
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    #m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    #den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = exp(-β*e[i]) - λ*exp(-β*e[j])
        #up = exp(-β*e[i]-m) - λ*exp(-β*e[j]-m)
        up = up * A[i,j] * B[j,i]
        down = z - e[j] + e[i]
        num += up/down
    end
    return num/den
end


"""
Masubara frequency G(iωn)/iωn = ∑_mn Cmn(iωn)/(Em - En): defalt type = :b
"""
function Masubara_freq_GFdivOmega(n::Integer, A::AbstractMatrix,B::AbstractMatrix,
                                ψ::CMPS, W::CMPO, β::Real)
    if n == 0 @error "Error: n should not be 0." end
    ωn = Masubara_freq(n,β,type=:b)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    #m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    #den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = A[i,j] * B[j,i] * diffaddexp(β,e[i],e[j])
        down = 1.0im * ωn - e[j] + e[i]
        num += up/down
    end
    return num/den
end


"""
    ∂ReG(iωn)/∂ωn
"""
function ∂ReG_∂ωn(n::Integer, A::AbstractMatrix,B::AbstractMatrix,
    ψ::CMPS, W::CMPO, β::Real)
    if n == 0 @error "ωn should not be 0." end
    λ = 1.0
    ωn = Masubara_freq(n,β,type=:b)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = exp(-β*e[i]) - λ*exp(-β*e[j]) 
        up = up * A[i,j] * B[j,i] * (e[i] - e[j]) * (-2ωn)
        down = ωn^2 + (e[i] - e[j])^2
        down = down^2
        num += up/down
    end
    return num/den
end

function ∂ReG_∂ωn_new(n::Integer, A::AbstractMatrix,B::AbstractMatrix,
    ψ::CMPS, W::CMPO, β::Real)
    if n == 0 @error "ωn should not be 0." end
    λ = 1.0
    ωn = Masubara_freq(n,β,type=:b)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = exp(-β*e[i]) - λ*exp(-β*e[j]) 
        up = up * A[i,j] * B[j,i] * (e[i] - e[j]) * (-2ωn)
        down = ωn^2 + (e[i] - e[j])^2
        down = down^2
        num += up/down

        up2 = ωn*exp(-β*e[j])
        down2 = ωn^2 + (e[i] - e[j])^2
        num += up2/down
    end
    return num/den
end



"""
spectral density: ρ(ω) = 2Imχ(ω) = -2ImG(ω)
"""
function spectral_density(ω::Real,A::AbstractMatrix,B::AbstractMatrix,
                             ψ::CMPS, W::CMPO, β::Real; η::Float64 = 0.001)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    #m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    #den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        res = exp(-β*e[i]) - exp(-β*e[j])
        #res = exp(-β*e[i]-m) - exp(-β*e[j]-m)
        res = res * A[i,j] * B[j,i] * delta(ω+e[i]-e[j],η)
        num += res
    end
    return 2π*num/den
end


"""
spectral density when A = B'
"""
function susceptibility(ω::Real,A::AbstractMatrix,
    ψ::CMPS, W::CMPO, β::Real; η::Float64 = 0.001)
    return spectral_density(ω, A, A',ψ, W, β; η=η)
end

function susceptibility(ω::Real, A::AbstractMatrix, B::AbstractMatrix,
    ψ::CMPS, W::CMPO, β::Real; η::Float64 = 0.001)
    if A' == B
        return spectral_density(ω, A, B, ψ, W, β; η=η)
    else
        @error "A' == B or use spectral_density function"
    end
end


"""
structure factor(spectral representation)
"""
function structure_factor(ω::Real, A::AbstractMatrix,B::AbstractMatrix,
                        ψ::CMPS, W::CMPO, β::Real; η::Float64 = 0.001)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    #m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    #den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        num += exp(-β*e[i])*A[i,j]*B[j,i]*delta(ω+e[i]-e[j], η)
        #num += exp(-β*e[i]-m)*A[i,j]*B[j,i]*delta(ω+e[i]-e[j], η)
    end
    return num/den * 2π
end

#end  # module PhysicalObservables
