#module PhysicalObservables

#using Reexport: @reexport

#include("SetupStruct.jl")
#@reexport using .SetupStruct

#export FreeEnergy
"""
Free energy
"""
function FreeEnergy(ψ::cmps, W::cmps, β::Real)
    Hψ = myprod(W,ψ)
    res = log(myinnerprod(ψ, Hψ ,β))- log(myinnerprod(ψ,ψ,β))
    return -1/β * res
end

function OptimFreeEnergy(x::Array{Float64,3})
    ψ = cmps(x[:,:,1], x[:,:,2])
    return FreeEnergy(ψ,W,β)
end
#end  # module PhysicalObservables
