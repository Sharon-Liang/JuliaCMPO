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
function TFIsing(J::Real, Γ::Real)
    return cmpo(Γ*pauli('x'), √J*pauli('z'), √J*pauli('z'), zeros(2,2))
end

"""
Free energy
"""
function FreeEnergy(ψ::cmps, W::cmpo, β::Real)
    Hψ = myprod(W,ψ)
    res = log(myinnerprod(ψ, Hψ ,β))- log(myinnerprod(ψ,ψ,β))
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
