"""
    PhysModel

Data structure of a physical model.
...
# Arguments 
- `cmpo`: the local transfer matrix.
- `phy_dim`: bond dimension of physical legs.
- `vir_dim`: virtual bond dimension of the horizontal legs.
- `Ut`: unitary transformation that transform a CMPO to its transpose. 
        Note that in general, `Ut` not necessarily be unitary, but be invertible, 
        however, in the limiting cases we know right now, they are unitary.
...
"""
@with_kw struct PhysModel{Tm<:AbstractCMPO, Tu}
    Tmatrix::Tm 
    phy_dim::Int64  
    vir_dim::Int64  
    Ut::Tu
end


"""
    generalUt(vD)

General form of `Ut` matrix ``Ut = [0 I; I 0]`` of a CMPO with virtual bond dimension `vD`.
"""
function generalUt(vD::Integer)
    zero_block = zeros(vD, vD)
    eye_block  = oneunit(zero_block)
    col1 = cat(zero_block, eye_block, dims = 1)
    col2 = cat(eye_block, zero_block, dims = 1)
    return cat(col1, col2, dims = 2), 2*vD
end

"""
    expand_cmpo(o::AbstractCMPO)

Expand a CMPO.
"""
function expand_cmpo(o::AbstractCMPO)
    R = √(0.5) * cat(o.R, o.L, dims = 3)
    L = √(0.5) * cat(o.L, o.R, dims = 3)

    if typeof(o.P) <: AbstractMatrix
        oP = reshape(o.P, size(o.P,1), size(o.P,2), 1, 1)
    else
        oP = o.P
    end

    pl = cat(oP, zeros(eltype(oP), size(oP)), dims = 3)
    pr = cat(zeros(eltype(oP), size(oP)), ein"ijkl->ijlk"(oP), dims = 3)
    P = cat(pl, pr, dims = 4)
    
    return CMPO(o.Q, R, L, P)
end