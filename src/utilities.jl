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
    diagm(v) function for CuVector
"""
LinearAlgebra.diagm(v::CuVector) = ein"i->ii"(v)
#end  # module utilities
