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

#end  # module PhysicalObservables
