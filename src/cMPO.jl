module cMPO

using Reexport: @reexport
using LinearAlgebra
using Zygote
using Optim
using Random; Random.seed!()

include("SetupStruct.jl")
@reexport using .SetupStruct

include("Models.jl")
include("PhysicalObservables.jl")

@show W |> typeof


end # module
