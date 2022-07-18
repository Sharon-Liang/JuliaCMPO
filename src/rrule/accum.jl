gentup(struct_T) = NamedTuple{( fieldnames(struct_T)...,), Tuple{(fieldtype(struct_T,i) for i=1:fieldcount(struct_T))...}}
@generated function to_named_tuple(x)
    nt = Expr(:quote, gentup(x))
    tup = Expr(:tuple)
    for i=1:fieldcount(x)
        push!(tup.args, :(getfield(x, $i)) )
    end
    return :($nt($tup))
end

function Zygote.accum(x::T,y::T) where T<:AbstractCMPS
    xtup, ytup = map(to_named_tuple, (x, y))
    res = Zygote.accum(xtup, ytup)
    return res
end

function Zygote.accum(x::AbstractCMPS{T,S,U}, y::NamedTuple{(:Q, :R), Tuple{S1,U1}}) where {T,S,U, S1,U1}
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