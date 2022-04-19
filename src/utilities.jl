#module utilities

"""
    Pauli matrices
"""
@enum PauliMatrixName PX PY iPY PZ PPlus PMinus
function pauli(T::DataType, name::PauliMatrixName)
    One = one(T)
    Zero = zero(T)
    if name == PX return [Zero One; One Zero]
    elseif name == iPY return [Zero One; -One Zero]
    elseif name == PZ return [One Zero; Zero -One]
    elseif name == PPlus return [Zero One; Zero Zero]
    elseif name == PMinus return [Zero Zero; One Zero]
    else return [Zero -One*im; One*im Zero] #iPY
    end
end
pauli(name::PauliMatrixName) = pauli(Float64, name)


"""
Gaussian approximate delta function
"""
function delta(x::Real, η::Real)
    num = η
    den = (x^2 + η^2) * π
    return num/den
end


"""
    Masubara frequency
"""
@enum OperatorType Bose Fermi
function Masubara_freq(n::Int64, β::Real; type::OperatorType = Bose)
    type == Bose ?  N = 2n : N = 2n + 1
    return N*π/β
end


"""
    Symmetrize a matrix `M`
"""
function symmetrize(M::AbstractMatrix)
   return (M + M')/2 #|> Hermitian
end


""" `symeigen` : manually symmetrize M before the eigen decomposition
    For CUDA dense matrix eigen solver `CUSOLVER.syevd!` and `CUSOLVER.heevd!`:
        'N'/'V': return eigenvalues/both eigenvalues and eigenvectors 
        'U'/'L': Upper/Lower triangle of `M` is stored.
"""
function symeigen(M::AbstractMatrix)
    e, v = M |> symmetrize |> eigen
    return e, v
end

for elty in (:Float32, :Float64)
    @eval begin
        function symeigen(M::CuMatrix{$elty})
            e, v = CUSOLVER.syevd!('V', 'U', symmetrize(M))
            return e, v
        end
    end
end

for elty in (:ComplexF32, :ComplexF64)
    @eval begin
        function symeigen(M::CuArray{$elty})
            e, v = CUSOLVER.heevd!('V', 'U', symmetrize(M))
            return e, v
        end
    end
end

"""
    `ln(Tr[exp(A)])` function
"""
function logtrexp(A::AbstractMatrix)
    val, _ = symeigen(A)
    return logsumexp(val)
end

"""
    `consist_diagm`: diagm(v) function that returns the same type of Matrix as input `v`
"""
consist_diagm(v::AbstractVector) = diagm(v)
consist_diagm(v::CuVector) = CuMatrix(CUDA.@allowscalar CUDA.diagm(v))
#end  # module utilities
