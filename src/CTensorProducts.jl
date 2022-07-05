"""Products of cmps and cmpo"""
function *(sl::AbstractCMPS{T, S, U}, sr::AbstractCMPS{T, S, U}) where {T, S, U}
    return CMPSMatrix(sl, sr)
end

function *(o::AbstractCMPO{T, S, U, V}, s::AbstractCMPS{T, S, U}) where {T, S, U, V}
    oi = convert(S, Matrix{T}(I,size(o.Q)))
    si = convert(S, Matrix{T}(I,size(s.Q)))
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return CMPS_generate(Q, R)
end

function *(s::AbstractCMPS{T, S, U}, o::AbstractCMPO{T, S, U, V}) where {T, S, U, V} 
    oi = convert(S, Matrix{T}(I,size(o.Q)))
    si = convert(S, Matrix{T}(I,size(s.Q)))
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return CMPS_generate(Q, R)
end

function *(ol::AbstractCMPO{T, S, U, V}, or::AbstractCMPO{T, S, U, V}) where {T, S, U, V}
    li = convert(S, Matrix{T}(I,size(ol.Q)))
    ri = convert(S, Matrix{T}(I,size(or.Q)))
    Q = ol.Q ⊗ ri + li ⊗ or.Q + ol.L ⊗ or.R
    L = li ⊗ or.L + ol.L ⊗ or.P
    R = ol.R ⊗ ri + ol.P ⊗ or.R
    P = ol.P ⊗ or.P
    return CMPO_generate(Q,R,L,P)
end