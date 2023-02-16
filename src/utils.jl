"""
    symmetrize(M)

Symmetrize `M`.
"""
symmetrize(M) = (M + M')/2


"""
    diagm(v::CuVector)

Construct a square matrix of type `CuMatrix` from `v` with its elements as diagonals. 
"""
LinearAlgebra.diagm(v::CuVector) = ein"i->ii"(v)


"""
    pauli([T=Float64,] name::PauliMatrixName)

Generate Pauli matrix `name`, with element type `T`.
"""
function pauli(::Type{T}, name::PauliMatrixName) where T
    One = one(T)
    Zero = zero(T)
    if name == PX 
        return [Zero One; One Zero]
    elseif name == iPY 
        return [Zero One; -One Zero]
    elseif name == PZ 
        return [One Zero; Zero -One]
    elseif name == PPlus 
        return [Zero One; Zero Zero]
    elseif name == PMinus 
        return [Zero Zero; One Zero]
    else #PY
        return [Zero -One*im; One*im Zero]
    end
end
pauli(name::PauliMatrixName) = pauli(Float64, name)



#=
### *Eigensolver* 
=#
"""
    eigensolver(A<:AbstractMatrix)

Computes the eigenvalue decomposition of a hermitian matrix `A`. The eigenvalues are sorted according to the smallest real part by defalt.

For a CUDA dense matrix, eigen solver `CUSOLVER.syevd!` and `CUSOLVER.heevd!` are used, where:

* 'N'/'V': return eigenvalues/both eigenvalues and eigenvectors.

* 'U'/'L': Upper/Lower triangle of `A` is stored. 

Note:
* `A` is explicitly symmetrized before calling `eigen` function.
* `A` itself is changed when `A` is a `CuMatrix`.
"""
function eigensolver end

eigensolver(A::AbstractArray; kwargs...) = eigen(symmetrize(A), kwargs...)

for elty in (:Float32, :Float64)
    @eval eigensolver(A::CuMatrix{$elty}) = CUSOLVER.syevd!('V', 'U', symmetrize(A))
end

for elty in (:ComplexF32, :ComplexF64)
    @eval eigensolver(A::CuMatrix{$elty}) = CUSOLVER.heevd!('V', 'U', symmetrize(A))
end


#=
### *logtrexp* 
=#
"""
    logtrexp(t, A)

Calculate ``lnTr(e^{tA})``, where ``t`` is a real number and ``A`` is a hermitian matrix.
"""
function logtrexp(t::Real, A::AbstractMatrix)
    vals, _ = eigensolver(A)
    return logsumexp(t*vals)
end



#=
### *Optim Functions* 
=#
#Reference: https://github.com/baggepinnen/FluxOptTools.jl/blob/master/src/FluxOptTools.jl
veclength(grads::Zygote.Grads) = sum(length(grads[p]) for p in grads.params)
veclength(params::Zygote.Params) = sum(length, params.params)
veclength(x) = length(x)

Base.zeros(grads::Zygote.Grads) = zeros(veclength(grads))
Base.zeros(pars::Zygote.Params) = zeros(veclength(pars))


"""
    optim_functions(loss, pars::Zygote.Params)

Generate vectorized parameter `p₀`, loss function and gradient function.
"""
function optim_functions(loss, pars::Zygote.Params)
    grads = Zygote.gradient(loss, pars)
    p₀ = zeros(pars)
    copy!(p₀, pars)

    # Gradient function
    g! = function (g,w)
        copy!(pars, w)
        grads = Zygote.gradient(loss, pars)
        copy!(g, grads)
    end

    # Preprocessed loss function
    f = function (w)
        copy!(pars, w)
        loss()
    end
    return p₀, f, g!
end


#=
### *Solver Functions* 
=#
"""
    solver_function(p::Processor)

Generate the solver function according to processor used.
"""
function solver_function(p::Processor)
    p == CPU ? solver = cpu_solver : solver = gpu_solver
    return solver
end


"""
    cpu_solver
"""
cpu_solver(a)= "Not Implemented" 
cpu_solver(a::Array) = a
cpu_solver(a::CuArray) = Array(a)

cpu_solver(a::CMPS{T,S,U}) where {T, S<:Matrix, U<:Array} = a

function cpu_solver(a::CMPS{T,S,U}) where {T, S<:CuMatrix, U<:CuArray} 
    @unpack Q, R = a
    Q = Matrix(Q)
    R = Array(R)
    return CMPS(Q, R)
end

cpu_solver(a::CMPO{T,S,U,V}) where {T, S<:Matrix, U<:Array, V<:Array} = a

function cpu_solver(a::CMPO{T,S,U,V}) where {T, S<:CuMatrix, U<:CuArray, V<:CuArray} 
    @unpack Q, R, L, P = a
    Q = Matrix(Q)
    R = Array(R)
    L = Array(L)
    P = Array(P)
    return CMPO(Q, R, L, P)
end


cpu_solver(expr::Function, M::Array...) = expr(M...)

function cpu_solver(expr::Function, M::CuArray...)
    M1 = Vector{Array}(undef, lastindex(M))
    for i in eachindex(M)
        M1[i] = convert(Array, M[i])
    end
    return expr(M1...)
end



"""
    gpu_solver
"""
gpu_solver(a)= "Not Implemented" 

gpu_solver(a::Array) = CuArray(a)
gpu_solver(a::CuArray) = a

gpu_solver(a::CMPS{T,S,U}) where {T, S<:CuMatrix, U<:CuArray} = a

function gpu_solver(a::CMPS{T,S,U}) where {T, S<:Matrix, U<:Array} 
    @unpack Q, R = a
    Q = CuMatrix(Q)
    R = CuArray(R)
    return CMPS(Q, R)
end

gpu_solver(a::CMPO{T,S,U,V}) where {T, S<:CuMatrix, U<:CuArray, V<:CuArray} = a

function gpu_solver(a::CMPO{T,S,U,V}) where {T, S<:Matrix, U<:Array, V<:Array} 
    @unpack Q, R, L, P = a
    Q = CuMatrix(Q)
    R = CuArray(R)
    L = CuArray(L)
    P = CuArray(P)
    return CMPO(Q, R, L, P)
end

gpu_solver(expr::Function, M::CuArray...) = expr(M...)

function gpu_solver(expr::Function, M::Array...)
    M1 = Vector{CuArray}(undef, lastindex(M))
    for i in eachindex(M)
        M1[i] = convert(CuArray, M[i])
    end
    return expr(M1...)
end