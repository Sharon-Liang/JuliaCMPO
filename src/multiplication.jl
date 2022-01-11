
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
function *(sl::CMPS, sr::CMPS)
    li = Matrix{eltype(sl.Q)}(I,size(sl.Q))
    ri = Matrix{eltype(sr.Q)}(I,size(sr.Q))
    K = li ⊗ sr.Q + sl.Q ⊗ ri + sl.R ⊗ sr.R
    return -K
end

function *(o::CMPO, s::CMPS) 
    oi = Matrix{eltype(o.Q)}(I,size(o.Q))
    si = Matrix{eltype(s.Q)}(I,size(s.Q))
    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
    R = o.R ⊗ si + o.P ⊗ s.R
    return CMPS(Q, R)
end

function *(s::CMPS, o::CMPO) 
    oi = Matrix{eltype(o.Q)}(I,size(o.Q))
    si = Matrix{eltype(s.Q)}(I,size(s.Q))
    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
    R = si ⊗ o.L + s.R ⊗ o.P
    return CMPS(Q, R)
end

function *(ol::CMPO, or::CMPO) 
    li = Matrix{eltype(ol.Q)}(I,size(ol.Q))
    ri = Matrix{eltype(or.Q)}(I,size(or.Q))
    Q = ol.Q ⊗ ri + li ⊗ or.Q + ol.L ⊗ or.R
    L = li ⊗ or.L + ol.L ⊗ or.P
    R = ol.R ⊗ ri + ol.P ⊗ or.R
    P = ol.P ⊗ or.P
    return CMPO(Q,R,L,P)
end