using Reexport: @reexport
using LinearAlgebra
using Zygote
using Optim
using Random; Random.seed!()

include("cMPO.jl")
@reexport using .cMPO
include("Models.jl")
include("PhysicalObservables.jl")

function randcmps(χ::Integer; D = 2, hermition = true)
    Q = rand(χ, χ)
    if D==2
        R = rand(χ, χ)
        if hermition
            Q = symmetrize(Q)
            R = symmetrize(R)
        end
    else
        R = rand(χ, χ, D-1)
    end

    return cmps(Q,R)
end

s = randcmps(2)
Hs = myprod(W, s)
TrExp(s.Q,1)
myinnerprod(Hs, Hs, 1)
