Zygote.@adjoint Array(x::CuArray) = Array(x), dy->(CuArray(dy),)

include("./logtrexp_Full_ED.jl")
include("./logtrexp_simple_FTLM.jl")
include("./logtrexp_orthogonalized_FTLM.jl")
include("./logtrexp_replaced_FTLM.jl")
include("./logtrexp_FullSampling_FTLM.jl")