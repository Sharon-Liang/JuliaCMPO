#Constructor adjoint
Zygote.@adjoint Array(x::CuArray) = Array(x), dy->(CuArray(dy),)

Zygote.@adjoint CMPS(Q, R) = CMPS(Q, R), ȳ->(ȳ.Q, ȳ.R)
Zygote.@adjoint CuCMPS(Q, R) = CuCMPS(Q, R), ȳ->(ȳ.Q, ȳ.R)

Zygote.@adjoint CMPSMatrix(ψl, ψr) = CMPSMatrix(ψl, ψr), ȳ->(ȳ.ψl, ȳ.ψr)
Zygote.@adjoint CTensor(x::T) where T<:CMPSMatrix = CTensor(x), ȳ->(CTensor(CMPSMatrix(ȳ.ψl, ȳ.ψr)), )
Zygote.@adjoint CuCTensor(x::T) where T<:CMPSMatrix = CuCTensor(x), ȳ->(CuCTensor(CMPSMatrix(ȳ.ψl, ȳ.ψr)),)

#accum function
gentup(struct_T) = NamedTuple{( fieldnames(struct_T)...,), Tuple{(fieldtype(struct_T,i) for i=1:fieldcount(struct_T))...}}
@generated function to_named_tuple(x)
    nt = Expr(:quote, gentup(x))
    tup = Expr(:tuple)
    for i=1:fieldcount(x)
        push!(tup.args, :(getfield(x, $i)) )
    end
    return :($nt($tup))
end

function Zygote.accum(x::AbstractCMPS,y::AbstractCMPS)
    xtup, ytup = map(to_named_tuple, (x, y))
    res = Zygote.accum(xtup, ytup)
    return res
end

function Zygote.accum(x::AbstractCMPS, y::NamedTuple{(:Q, :R), Tuple{Tq,Tr}}) where {Tq, Tr}
    xtup = to_named_tuple(x)
    res = Zygote.accum(xtup, y)
    return res
end

function Zygote.accum(x::T, y::T, zs::T...) where T<:AbstractCMPS
    xtup, ytup = map(to_named_tuple, (x, y))
    zstup = map(to_named_tuple, zs)
    res = Zygote.accum(xtup, ytup,zstup...)
    return res
end

function _sum_expr_w1_w2(expr::Function, vals, w1, w2)
    Λ = map(expr, vals)
    return reduce(+, map(*, Λ, w1, w2))
end

include("./logtrexp_full_ed.jl")
include("./logtrexp_standard_ftlm.jl")
include("./logtrexp_replaced_ftlm.jl")
include("./logtrexp_orthogonalized_ftlm.jl")