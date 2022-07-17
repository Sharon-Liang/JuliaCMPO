Zygote.@adjoint Array(x::CuArray) = Array(x), dy->(CuArray(dy),)

function ChainRules.rrule(::Type{CMPSMatrix}, ψl::AbstractCMPS{T, S, U}, 
                        ψr::AbstractCMPS{T, S, U}) where {T,S,U}
    CMPSMatrix_pullback(ȳ) = ChainRules.NoTangent(), ȳ.ψl, ȳ.ψr
    return CMPSMatrix(ψl, ψr), CMPSMatrix_pullback
end

Zygote.@adjoint CTensor(x::T) where T<:CMPSMatrix = CTensor(x), ȳ->(CTensor(CMPSMatrix(ȳ.ψl, ȳ.ψr)), )
Zygote.@adjoint CuCTensor(x::T) where T<:CMPSMatrix = CuCTensor(x), ȳ->(CuCTensor(CMPSMatrix(ȳ.ψl, ȳ.ψr)),)

include("./logtrexp_Full_ED.jl")
include("./logtrexp_simple_FTLM.jl")
include("./logtrexp_orthogonalized_FTLM.jl")
include("./logtrexp_replaced_FTLM.jl")
include("./logtrexp_FullSampling_FTLM.jl")
