#module PhysicalObservables
#include("Setup.jl")

function pauli(char::Char)
    if char=='x' return [0. 1.; 1. 0.]
    elseif char=='y' return [0. -1im; 1im 0.]
    elseif char=='z' return [1. 0.; 0. -1.]
    elseif char=='+' return [0. 1.; 0. 0.]
    elseif char=='-' return [0. 0.; 1. 0.]
    else
        return printf("Warnning: the input should be 'x','y','z','+','_'.")
    end
end

"""
NN Transvers field Ising model
    H = ∑ J Zi Zj + ∑ Γ Xi
"""
function TFIsing(J::Real, Γ::Real; field = 0.0)
    h = field * pauli('z')
    return cmpo(Γ*pauli('x') + h, √J*pauli('z'), √J*pauli('z'), zeros(2,2))
end

"""
Free energy
"""
function FreeEnergy(ψ::cmps, W::cmpo, β::Real)
    Hψ = W * ψ
    res = LogTrExp(ψ * Hψ, β)- LogTrExp(ψ * ψ, β)
    return -1/β * res
end

function OptimFreeEnergy(x::Array{Float64,3}, W::cmpo, β::Real)
    ψ = cmps(x[:,:,1], x[:,:,2])
    return FreeEnergy(ψ,W,β)
end

function OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W::cmpo, β::Real)
    ψ = cmps(x[:,:,1], x[:,:,2])
    grad = gradient(ψ -> FreeEnergy(ψ, W, β), ψ)[1]
    (r,c) = size(grad.Q)
    for i = 1:r, j = 1:c
        gx[i,j,1] = grad.Q[i,j]
        gx[i,j,2] = grad.R[i,j]
    end
end

"""
The thermal average of local opeartors
"""
function Thermal_average(ψ::cmps, W::cmpo, Op::AbstractArray, β::Real)
    eye = Matrix(1.0I, size(ψ.Q))
    Op = kron(eye, kron(Op, eye))
    K = ψ * (W * ψ)
    K = symmetrize(K)
    vals, U = eigen(K)
    m = maximum(vals)
    Op = U' * Op * U
    den = exp.(β* (vals .- m)) |> sum
    num = exp.(β * (vals .- m)) .* diag(Op) |> sum
    return num/den
end

"""
The local two-time correlation functions
"""
function Correlation_2time(A::AbstractArray,B::AbstractArray,ψ::cmps, W::cmpo, β::Real, τ::Number)
       eye = Matrix(1.0I, size(ψ.Q))
       A = kron(eye, kron(A, eye))
       B = kron(eye, kron(B, eye))
       K = ψ * (W * ψ)
       K = symmetrize(K)
       vals, U = eigen(K)
       m = maximum(vals)
       A = U' * A * U |> symmetrize
       B = U' * B * U |> symmetrize
       den = exp.(β* (vals .- m)) |> sum
       num = 0.0
       for i = 1: length(vals), j = 1: length(vals)
           num += exp(β*(vals[i]-m)-τ*(vals[i] - vals[j])) * A[i,j] * B[j,i]
       end
       return num/den
end


function OptimDiff(x::Array{Float64,3})
    ψ = cMPS(x[:,:,1], x[:,:,2])
    return difference(ψ,Hψ)
end

function OptimDiff!(gx::Array{Float64,3}, x::Array{Float64,3})
    ψ = cMPS(x[:,:,1], x[:,:,2])
    grad = gradient(ψ -> difference(ψ,Hψ), ψ)[1]
    (r,c) = size(grad.Q)
    for i = 1:r, j = 1:c
        gx[i,j,1] = grad.Q[i,j]
        gx[i,j,2] = grad.R[i,j]
    end
end
#end  # module PhysicalObservables
