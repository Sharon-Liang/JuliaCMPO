import Base:*
"""
    *(a::AbstractCTensor, b::AbstractCTensor)
"""
*(sl::AbstractCMPS, sr::AbstractCMPS) = CMPSMatrix(sl, sr)


function *(o::AbstractCMPO, s::AbstractCMPS)
    oi = Matrix{eltype(o.Q)}(I, size(o.Q)) 
    si = Matrix{eltype(s.Q)}(I, size(s.Q))
    oi = convert(typeof(o.Q), oi)
    si = convert(typeof(s.Q), si)
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return cmps_generate(Q, R)
end

function *(s::AbstractCMPS, o::AbstractCMPO)
    oi = Matrix{eltype(o.Q)}(I, size(o.Q))
    si = Matrix{eltype(s.Q)}(I, size(s.Q))
    oi = convert(typeof(o.Q), oi)
    si = convert(typeof(s.Q), si)
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return cmps_generate(Q, R)
end

function *(ol::AbstractCMPO, or::AbstractCMPO)
    li, ri = oneunit(ol.Q), oneunit(or.Q)
    Q = ol.Q ⊗ ri + li ⊗ or.Q + ol.L ⊗ or.R
    L = li ⊗ or.L + ol.L ⊗ or.P
    R = ol.R ⊗ ri + ol.P ⊗ or.R
    P = ol.P ⊗ or.P
    return cmpo_generate(Q,R,L,P)
end