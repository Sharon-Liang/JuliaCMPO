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
    res = logtrexp(-β*(ψ*W*ψ))- logtrexp(-β*(ψ * ψ))
    return -1/β * res
end

function free_energy(param::Array{Float64,3}, W::cmpo, β::Real)
    ψ = cmps(param[:,:,1], param[:,:,2])
    free_energy(ψ, W, β)
end


"""
The thermal average of local opeartors
"""
function thermal_average(ψ::cmps, W::cmpo, Op::AbstractArray, β::Real)
    eye = Matrix(1.0I, size(ψ.Q))
    Op = eye ⊗ Op ⊗ eye
    K = ψ * W * ψ |> symmetrize
    vals, U = eigen(K)
    m = maximum(vals)
    Op = U' * Op * U
    den = exp.(β * (vals .- m)) |> sum
    num = exp.(β * (vals .- m)) .* diag(Op) |> sum
    return num/den
end

"""
The local two-time correlation functions
"""
function correlation_2time(A::AbstractArray,B::AbstractArray,ψ::cmps, W::cmpo, β::Real, τ::Number)
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

#end  # module PhysicalObservables
