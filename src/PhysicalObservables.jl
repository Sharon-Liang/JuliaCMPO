#module PhysicalObservables
#include("Setup.jl")

"""
NN Transvers field Ising model
    H = ∑ J Zi Zj + ∑ Γ Xi
"""
function TFIsing(J::Real, Γ::Real; field='N')
    if field == 'N'
        h = zeros(2,2)
    else
        η = 1.e-5
        h = η .* pauli(field)
    end
    return cmpo(Γ*pauli('x')+h, √J*pauli('z'), √J*pauli('z'), zeros(2,2))
end

"""
Free energy
"""
function free_energy(ψ::cmps, W::cmpo, β::Real)
    K = ψ * W * ψ ; H = ψ * ψ
    res = logtrexp(-β*K)- logtrexp(-β*H)
    return -1/β * res
end

function free_energy(param::Array{Float64,3}, W::cmpo, β::Real)
    ψ = cmps(param[:,:,1], param[:,:,2])
    free_energy(ψ, W, β)
end


"""
The thermal average of local opeartors
"""
function thermal_average(Op::AbstractArray, ψ::cmps, W::cmpo, β::Real)
    eye = Matrix(1.0I, size(ψ.Q))
    Op = eye ⊗ Op ⊗ eye
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(-β*K)
    m = maximum(e)
    Op = v' * Op * v
    den = exp.(e .- m) |> sum
    num = exp.(e .- m) .* diag(Op) |> sum
    return num/den
end

"""
The local two-time correlation functions
"""
function correlation_2time(τ::Number, A::AbstractArray,B::AbstractArray,
                           ψ::cmps, W::cmpo, β::Real)
    eye = Matrix(1.0I, size(ψ.Q))
    A = eye ⊗ A ⊗ eye
    B = eye ⊗ B ⊗ eye
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        num += exp(-β*e[i]- m + τ*(e[i] - e[j])) * A[i,j] * B[j,i]
    end
    return num/den
end

function susceptibility(n::Integer, A::AbstractArray,B::AbstractArray,
                        ψ::cmps, W::cmpo, β::Real)
    # i ωn
    ωn = 2π * n/β  #bosion
    eye = Matrix(1.0I, size(ψ.Q))
    A = eye ⊗ A ⊗ eye
    B = eye ⊗ B ⊗ eye
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = exp(- β*e[j]-m) - exp(-β * e[i]-m)
        up = up * A[i,j] * B[j,i]
        down = 1im*ωn + e[i] - e[j]
        num += up/down
    end
    return num/den |> real
end


function imag_susceptibility(ω::Real,A::AbstractArray,B::AbstractArray,
                             ψ::cmps, W::cmpo, β::Real; η::Float64 = 0.05)
    eye = Matrix(1.0I, size(ψ.Q))
    A = eye ⊗ A ⊗ eye
    B = eye ⊗ B ⊗ eye
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = exp(-β*e[j]-m) - exp(-β*e[i]-m)
        up = up * A[i,j] * B[j,i] * η
        down = (ω + e[i] - e[j])^2 + η^2
        num += up/down
    end
    return -num/den
end

function NMR_relaxation(A::AbstractArray,B::AbstractArray,
                        ψ::cmps, W::cmpo, β::Real; η::Float64 = 0.05)
    eye = Matrix(1.0I, size(ψ.Q))
    A = eye ⊗ A ⊗ eye
    B = eye ⊗ B ⊗ eye
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = exp(-β*e[j]-m) - exp(-β*e[i]-m)
        up = up * A[i,j] * B[j,i] * η *(e[i]-e[j])
        down = (e[i] - e[j])^2 + η^2
        down = down^2
        num += up/down
    end
    return num/den * 4/β
end
#end  # module PhysicalObservables
