
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
