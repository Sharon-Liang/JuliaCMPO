#=
### *Enumerates*
=#

"""
    OperatorType Bose Fermi

Operator type of physical observables
"""
@enum OperatorType Bose Fermi


"""
    PauliMatrixName PX PY iPY PZ PPlus PMinus

Name of Pauli matrices
"""
@enum PauliMatrixName PX PY iPY PZ PPlus PMinus


#=
### *Struct* : *Tensors*
=#
""" 
    CMPS{T<:Number, S<:AbstractMatrix, U<:AbstractArray}

Data structure of continuous matrix product states (cMPS) local tensor. The structure of a cMPS local tensor is: 
```
--       --
| I + ϵQ  |
|         |
|    |    |
|   √ϵR   |
|    |    |
--       --
```
where `Q` is a `χ×χ` dimensional matrix, `R` is a `χ×χ×D` dimensional array, and `χ` is the bond dimension of the cMPS, `(D+1)` is the virtual bond dimension in spatial direction.
"""
@with_kw struct CMPS{T<:Number, S<:AbstractMatrix, U<:AbstractArray}
    Q::S 
    R::U 
    CMPS{T,S,U}(Q::AbstractMatrix{T}, R::AbstractArray{T}) where {T,S,U} = new(Q,R)
end
CMPS(Q::AbstractMatrix{T}, R::AbstractArray{T}) where {T} = CMPS{T,typeof(Q), typeof(R)}(Q,R)


"""
    CMPO{T<:Number, S<:AbstractMatrix, U<:AbstractArray, V<:AbstractArray}

Data structure of continuous matrix product operator (cMPO) local tensor. The structure of a cMPO local tensor is: 
```
--                 --
| I + ϵQ  -- √ϵL -- |
|                   |
|   |               |
| √ϵR        P      |
|   |               |
--                 --
```
where `Q` is a `d×d` dimensional matrix, `R` and `L` are `d×d×D` dimensional arrays, and `d` is the physical bond dimension, `(D+1)` is the virtual bond dimension in spatial direction.
"""
@with_kw struct CMPO{T<:Number, S<:AbstractMatrix, U<:AbstractArray, V<:AbstractArray} 
    Q::S 
    R::U 
    L::U 
    P::V  
    CMPO{T,S,U,V}(Q::AbstractMatrix{T}, R::AbstractArray{T}, L::AbstractArray{T}, P::AbstractArray{T}) where {T,S,U,V} = new(Q,R,L,P)
end
CMPO(Q::AbstractMatrix{T}, R::AbstractArray{T}, L::AbstractArray{T}, P::Array{T}) where {T} = CMPO{T, typeof(Q), typeof(R), typeof(P)}(Q,R,L,P)


#=
### *Struct* : *Options*
=#

#TODO