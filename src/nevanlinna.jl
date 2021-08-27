# N_freq
# ωn = [Masubara_freq(n, β) for n = 1: N_freq]
# Gn = [Masubara_GF(n,z,z, ψ, w, β) for n = 1: N_freq]
# θn = [invMobiusTransform(-g) for g in Gn]
# check the existence of the Nevanlinna interpolant of θn

function MobiusTransform(z::ComplexF32)
    return (z - 1.0im)/(z + 1.0im)
end

function invMobiusTransform(z::ComplexF32)
    return 1.0im*(1 + z)/(1 - z)
end

function pick_matrix(freq::Vector, val::Vector)
    dim = length(freq)
    if length(val) != dim
        @error DimensionMismatch("dimension of val must match freq")
    end
    pmat = zeros(ComplexF32, dim, dim)
    for i = 1:dim, j = 1:dim
        num = 1 - val[i]*val[j]'
        den = 1 - invMobiusTransform(freq[i])*invMobiusTransform(freq[j])'
        pmat[i,j] = num/den
    end
    return pmat
end

function isNevanlinnasolvable(freq::Vector, val::Vector)
    # check if the pick_matrix is positive semidefinite
    pmat = pick_matrix(freq, val)
    (r,c) = size(r,c)
    pmat += Matrix(1.0I, r, c) * 1.e-250
    return pmat |> isposdef
end

function isNevanlinnaunique(freq::Vector, val::Vector)
    return pick_matrix(freq, val) |> det |> iszero
end

