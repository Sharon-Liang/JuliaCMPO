include("setup_struct.jl")

"""
Free energy
"""
function F(ψ::cMPS, W::cMPO, β::Real)
    Hψ = myprod(W,ψ)
    res = log(myinnerprod(ψ, Hψ ,β))- log(myinnerprod(ψ,ψ,β))
    return -1/β * res
end
