#module setup
"""
    Flattern cMPS: cMPS <-> Array
"""
function toarray(ψ::CMPS)
    sq = size(ψ.Q)
    sr = size(ψ.R)  #sr: dimension of ψ.R array
    if length(sr) == 2  #sr = 2, ψ.R is a matrix
        Q = reshape(ψ.Q, sq[1],sq[2],1)
        R = reshape(ψ.R, sr[1],sr[2],1)
    elseif length(sr) > 2
        Q = reshape(ψ.Q, sq[1],sq[2],1)
        R = ψ.R
    else
        @error "Improper CMPS"
    end
    return cat(Q,R,dims=3)
end

function tovector(ψ::CMPS)
    arr = ψ |> toarray
    dim = size(arr)
    return vec(arr), dim
end

function tocmps(A::AbstractArray{T,3} where T)
    d = size(A)[3]
    if d == 2
        return CMPS(A[:,:,1],A[:,:,2])
    else
        return CMPS(A[:,:,1],A[:,:,2:end])
    end
end

function tocmps(V::Vector{T} where T, dim::Tuple)
    arr = reshape(V, dim)
    return tocmps(arr)
end



"""
    Normalize a cMPS, i.e. ⟨ψ|ψ⟩ = 1
"""
#function normalize(s::CMPS, β::Real)
#    T = -β*(s*s) |> trexp
#    λ = (T.max + log(T.res))/β
#    λ = λ/2
#    eye = Matrix(λ*I,size(s.Q))
#    Q = s.Q - eye
#    return CMPS(Q, s.R)
#end


"""
    Logarithm of the overlap of two cMPS:
        log_overlap = ln(⟨ψl|ψr⟩)
"""
function log_overlap(sl::CMPS, sr::CMPS, β::Real)
    K = sl * sr
    return logtrexp(-β * K)
end


"""
    transpose of a CMPO
"""
function transpose(o::CMPO)
    Q = o.Q
    R = o.L
    L = o.R

    if length(size(o.P)) == 2
        oP = reshape(o.P, size(o.P)[1], size(o.P)[2], 1, 1)
    else
        oP = o.P
    end
    P = ein"ijkl->ijlk"(oP)

    return CMPO(Q,R,L,P)
end


"""
    Hermitian conjuate (adjoint) of a CMPO
"""
function adjoint(o::CMPO)
    Q = conj.(o.Q)
    R = conj.(o.L)
    L = conj.(o.R)
    if length(size(o.P)) == 2
        oP = reshape(o.P, size(o.P)[1], size(o.P)[2], 1, 1)
    else
        oP = o.P
    end
    P = conj.(ein"ijkl->ijlk"(oP))
    return CMPO(Q,R,L,P)
end


"""
    isequal(a::CMPO, b::CMPO)
"""
function isequal(a::CMPO, b::CMPO)
    res = Vector{Bool}(undef, 4)
    res[1] = all(isequal.(a.Q, b.Q))
    res[2] = all(isequal.(a.R, b.R))
    res[3] = all(isequal.(a.L, b.L))
    res[4] = all(isequal.(a.P, b.P))
    return all(res)
end


"""
    Perform a unitary transformation in the imaginary-time direction
    (i.e. One physical bond of cMPO and virtual bond of cMPS)
    if U is a square matrix, this is a guage transformation
"""
function project(s::CMPS, u::AbstractMatrix)
    Q = u' * s.Q * u
    if length(size(s.R)) == 2
        sR = reshape(s.R, size(s.R)[1], size(s.R)[2], 1)
    else
        sR = s.R
    end
    R = ein"(ip,pql),qj -> ijl"(u', sR, u)
    return CMPS(Q, R)
end


function project(o::CMPO, u::AbstractMatrix)
    Q = u' * o.Q * u
    if length(size(o.R)) == 2
        oR = reshape(o.R, size(o.R)[1], size(o.R)[2], 1)
        oL = reshape(o.L, size(o.L)[1], size(o.L)[2], 1)
        oP = reshape(o.P, size(o.P)[1], size(o.P)[2], 1, 1)
    else
        oR = o.R
        oL = o.L
        oP = o.P
    end
    R = ein"(ip,pql),qj -> ijl"(u', oR, u)
    L = ein"(ip,pql),qj -> ijl"(u', oL, u)
    P = ein"(ip,pqlm),qj -> ijlm"(u', oP, u)
    return CMPO(Q,R,L,P)
end


""" multiply a matrix to the left/right of the cMPS/CMPO
    --        --   --            --
    | 1 0 ... 0|   | I + dtau Q   |
    | 0        |   |              |
    | :        |   |       |      |
    | :    W   |   | sqrt(dtau) R |
    | 0        |   |       |      |
    --        --   --            --
"""
function *(Ut::AbstractMatrix, s::CMPS)
    if length(size(s.R)) == 2
        sR = reshape(s.R, size(s.R)[1], size(s.R)[2], 1)
    else
        sR = s.R
    end
    R = ein"mn, ijn -> ijm"(Ut, sR)
    return CMPS(s.Q, R)
end

function *(Ut::AbstractMatrix, o::CMPO)
    if length(size(o.R)) == 2
        oR = reshape(o.R, size(o.R)[1], size(o.R)[2], 1)
        oP = reshape(o.P, size(o.P)[1], size(o.P)[2], 1, 1)
    else
        oR = o.R
        oP = o.P
    end
    R = ein"mk, ijk -> ijm"(Ut, oR)
    P = ein"mk, ijkn -> ijmn"(Ut, oP) 
    return CMPO(o.Q, R, o.L, P)
end

function *(o::CMPO, Ut::AbstractMatrix)
    if length(size(o.L)) == 2
        oL = reshape(o.L, size(o.L)[1], size(o.L)[2], 1)
        oP = reshape(o.P, size(o.P)[1], size(o.P)[2], 1, 1)
    else
        oL = o.L
        oP = o.P
    end
    L = ein"ijk, km -> ijm"(oL, Ut)
    P = ein"ijnk, km -> ijnm"(oP, Ut) 
    return CMPO(o.Q, o.R, L, P)
end


""" 
    transform the cMPS to the gauge where Q is a diagonalized matrix
"""
function diagQ(s::CMPS; dim::Integer = 999)
    χ = min(size(s.Q)[1], dim)
    _, u = symeigen(s.Q)
    return project(s, u[:,end+1-χ:end])
end


#end module setup
