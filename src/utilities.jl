#module utilities
import Base: kron

"""
    Pauli matrices
"""
function pauli(T::DataType, symbol::Symbol)
    One = one(T)
    Zero = zero(T)
    if symbol==:x return [Zero One; One Zero]
    elseif symbol==:y return [Zero -One*im; One*im Zero]
    elseif symbol==:iy return [Zero One; -One Zero]
    elseif symbol==:z return [One Zero; Zero -One]
    elseif symbol==:+ return [Zero One; Zero Zero]
    elseif symbol==:- return [Zero Zero; One Zero]
    else
        error("The input should be :x,:y,:z,:+,:-, :iy.")
    end
end

function pauli(symbol::Symbol)
    return pauli(Float64, symbol)
end

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
function Masubara_freq(n::Int64, β::Real; type::Symbol=:b)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end


"""
    Symmetrize a matrix M
"""
function symmetrize(M::AbstractMatrix)
   return (M + M')/2 #|> Hermitian
end


""" symeigen
    manually symmetrize M before the eigen decomposition
"""
function symeigen(M::AbstractMatrix)
    return M |> symmetrize |> eigen
end


"""
    ln(Tr[exp(A)]) function
"""
function logtrexp(A::AbstractMatrix)
    val, _ = symeigen(A)
    return logsumexp(val)
end

#end  # module utilities
