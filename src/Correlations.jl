#module Correlations
"""
    The local two-time correlation functions: C(τ) = -G(τ) = ⟨A(τ)B(0)⟩
"""
function correlation_2time(τ::Number, A::AbstractMatrix,B::AbstractMatrix,
                            ψl::CMPS, ψr::CMPS, W::CMPO, β::Real)
    K = ψl * W * ψr 
    e, v = eigensolver(K)
    min = minimum(e); e = e .- min
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        num += exp(-β*e[i] + τ*(e[i] - e[j])) * A[i,j] * B[j,i]
    end
    return num/den
end
correlation_2time(τ::Number, A::AbstractMatrix, B::AbstractMatrix, ψ::CMPS, W::CMPO, β::Real) = 
    correlation_2time(τ, A, B, ψ, ψ, W, β)














"""
A useful function: f = e^(-b*e1) - e^(-b*e2) / (e2 - e1)
"""
function diffaddexp(b::Real, e1::Real, e2::Real)
    if abs(e2 - e1) < 1.e-10
        return exp(-b*e1) * b
    else
        num = exp(-b*e1) - exp(-b*e2)
        den = e2 - e1
        return num/den
    end
end

"""
Masubara frequency Green's functions: defalt type = :b
"""
function Masubara_freq_GF(n::Integer, A::AbstractMatrix,B::AbstractMatrix,
                        ψ::CMPS, W::CMPO, β::Real)
    λ = 1.0
    ωn = Masubara_freq(n,β,type=:b)
    K = ψ * W * ψ 
    e, v = eigensolver(K)
    min = minimum(e); e = e .- min
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    num = 0.0
    # sum i != j , add i = j term when ωn = 0
    if ωn != 0
        for i = 1: length(e), j = 1: length(e)
            up = exp(-β*e[i]) - λ*exp(-β*e[j])
            up = up * A[i,j] * B[j,i]
            down = 1.0im * ωn - e[j] + e[i]
            num += up/down
        end
    else
        for i = 1: length(e), j = 1: length(e)
            num -= A[i,j] * B[j,i]* diffaddexp(β,e[i],e[j])
        end
    end
    return num/den
end

#add spectrum Function

#end  # module Correlations
