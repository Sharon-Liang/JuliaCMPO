
"""
    multiplications of Arrays in CMPS and CMPO struct
"""
function ⊗(A::AbstractMatrix, B::AbstractMatrix)
    (r1, c1) = size(A)
    (r2, c2) = size(B)
    return reshape(ein"ij,kl->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::AbstractMatrix, B::AbstractArray{T,3} where T) 
    (r1, c1) = size(A)
    (r2, c2, d) = size(B)
    return reshape(ein"ij,klm->kiljm"(A, B), r1*r2, c1*c2, d)
end

function ⊗(A::AbstractArray{T,3} where T, B::AbstractMatrix)
    (r1, c1, d) = size(A)
    (r2, c2) = size(B)
    return reshape(ein"ijm,kl->kiljm"(A, B), r1*r2, c1*c2, d)
end

function ⊗(A::AbstractArray{T,3} where T, B::AbstractArray{T,3} where T)
    (r1, c1, d1) = size(A)
    (r2, c2, d2) = size(B)
    if d1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijm,klm->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::AbstractArray{T,4} where T, B::AbstractArray{T,3} where T)
    (r1, c1, d1, f) = size(A)
    (r2, c2, d2) = size(B)
    if f != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijnm,klm->kiljn"(A, B), r1*r2, c1*c2, d1)    
end

function ⊗(A::AbstractArray{T,3} where T, B::AbstractArray{T,4} where T)
    (r1, c1, d1) = size(A)
    (r2, c2, d2, f) = size(B)
    if d1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijm,klmn->kiljn"(A, B), r1*r2, c1*c2, f)    
end

function ⊗(A::AbstractArray{T,4} where T, B::AbstractArray{T,4} where T)
    (r1,c1,d1,f1) = size(A)
    (r2,c2,d2,f2) = size(B)
    if f1 != d2 @error "Dimension mismatch!" end
    return reshape(ein"ijpm,klmq->kiljpq"(A, B), r1*r2, c1*c2, d1, f2)    
end


"""multiplications of cmps and cmpo"""
function *(sl::AbstractCMPS{T}, sr::AbstractCMPS{T}) where T
    li = convert(typeof(sl.Q), Matrix{T}(I,size(sl.Q)))
    ri = convert(typeof(sr.Q), Matrix{T}(I,size(sr.Q)))
    K = li ⊗ sr.Q + sl.Q ⊗ ri + sl.R ⊗ sr.R
    return -K
end

function *(o::AbstractCMPO{T}, s::AbstractCMPS{T}) where T
    oi = convert(typeof(o.Q), Matrix{T}(I,size(o.Q)))
    si = convert(typeof(s.Q), Matrix{T}(I,size(s.Q)))
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return CMPS_generate(Q, R)
end

function *(s::AbstractCMPS{T}, o::AbstractCMPO{T}) where T 
    oi = convert(typeof(o.Q), Matrix{T}(I,size(o.Q)))
    si = convert(typeof(s.Q), Matrix{T}(I,size(s.Q)))
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return CMPS_generate(Q, R)
end

function *(ol::AbstractCMPO{T}, or::AbstractCMPO{T}) where T
    li = convert(typeof(ol.Q), Matrix{T}(I,size(ol.Q)))
    ri = convert(typeof(or.Q), Matrix{T}(I,size(or.Q)))
    Q = ol.Q ⊗ ri + li ⊗ or.Q + ol.L ⊗ or.R
    L = li ⊗ or.L + ol.L ⊗ or.P
    R = ol.R ⊗ ri + ol.P ⊗ or.R
    P = ol.P ⊗ or.P
    return CMPO_generate(Q,R,L,P)
end