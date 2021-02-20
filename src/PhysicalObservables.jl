module PhysicalObservables

using LinearAlgebra

include("ToolFunctions,jl")
include("SetupStruct.jl")

"""
Free energy
"""
function FreeEnergy(ψ::cMPS, W::cMPO, β::Real)
    Hψ = myprod(W,ψ)
    res = log(myinnerprod(ψ, Hψ ,β))- log(myinnerprod(ψ,ψ,β))
    return -1/β * res
end

end  # module PhysicalObservables
