Zygote.@adjoint Array(x::CuArray) = Array(x), dy->(CuArray(dy),)

Zygote.@adjoint CMPS(Q, R) = CMPS(Q, R), ȳ->(ȳ.Q, ȳ.R)
Zygote.@adjoint CuCMPS(Q, R) = CuCMPS(Q, R), ȳ->(ȳ.Q, ȳ.R)

Zygote.@adjoint CMPSMatrix(ψl, ψr) = CMPSMatrix(ψl, ψr), ȳ->(ȳ.ψl, ȳ.ψr)
Zygote.@adjoint CTensor(x::T) where T<:CMPSMatrix = CTensor(x), ȳ->(CTensor(CMPSMatrix(ȳ.ψl, ȳ.ψr)), )
Zygote.@adjoint CuCTensor(x::T) where T<:CMPSMatrix = CuCTensor(x), ȳ->(CuCTensor(CMPSMatrix(ȳ.ψl, ȳ.ψr)),)