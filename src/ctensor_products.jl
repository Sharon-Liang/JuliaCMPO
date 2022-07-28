import Base:*
"""
    *(a::AbstractCTensor, b::AbstractCTensor)
"""
*(sl::AbstractCMPS, sr::AbstractCMPS) = CMPSMatrix(sl, sr)


function *(o::AbstractCMPO, s::AbstractCMPS)
    oi, si = oneunit(o.Q), oneunit(s.Q)
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return CMPS_generate(Q, R)
end

function *(s::AbstractCMPS, o::AbstractCMPO)
    oi, si = oneunit(o.Q), oneunit(s.Q)
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return CMPS_generate(Q, R)
end

function *(ol::AbstractCMPO, or::AbstractCMPO)
    li, ri = oneunit(ol.Q), oneunit(or.Q)
    Q = ol.Q ⊗ ri + li ⊗ or.Q + ol.L ⊗ or.R
    L = li ⊗ or.L + ol.L ⊗ or.P
    R = ol.R ⊗ ri + ol.P ⊗ or.R
    P = ol.P ⊗ or.P
    return CMPO_generate(Q,R,L,P)
end