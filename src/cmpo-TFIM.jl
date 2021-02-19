using LinearAlgebra
using Zygote
using Optim
using Random; Random.seed!()

include("setup_struct.jl")
include("toolfunctions.jl")
include("models.jl")
include("phyobservables.jl")

β = 1.
