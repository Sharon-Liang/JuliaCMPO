
"""
    ⊗(A, B)

multiplications of Arrays in `AbstractCMPS` and `AbstractCMPO` struct
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
    r1, c1, _ = size(A)
    r2, c2, _ = size(B)
    return reshape(ein"ijm,klm->kilj"(A, B), r1*r2, c1*c2)
end

function ⊗(A::AbstractArray{T,4} where T, B::AbstractArray{T,3} where T)
    r1, c1, d1, _ = size(A)
    r2, c2, _ = size(B)
    return reshape(ein"ijnm,klm->kiljn"(A, B), r1*r2, c1*c2, d1)    
end

function ⊗(A::AbstractArray{T,3} where T, B::AbstractArray{T,4} where T)
    r1, c1, _ = size(A)
    r2, c2, _, f = size(B)
    return reshape(ein"ijm,klmn->kiljn"(A, B), r1*r2, c1*c2, f)    
end

function ⊗(A::AbstractArray{T,4} where T, B::AbstractArray{T,4} where T)
    r1,c1,d1,_ = size(A)
    r2,c2,_,f2 = size(B)
    return reshape(ein"ijpm,klmq->kiljpq"(A, B), r1*r2, c1*c2, d1, f2)    
end
