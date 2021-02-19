struct cMPS{T<:AbstractArray}
    Q::T
    R::T
end

struct cMPO{T<:AbstractArray}
    Q::T  # onsite
    R::T  # interaction
    L::T  # interaction
    P::T  # long-range
end

function toarray(ψ::cMPS)
    # size(Q) == size(R)
    (r,c) = size(ψ.Q)
    x = zeros(r,c,2)
    x[:,:,1] = ψ.Q
    x[:,:,2] = ψ.R
    return x
end
