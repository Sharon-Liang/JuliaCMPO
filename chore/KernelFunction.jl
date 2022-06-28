using ChainRulesCore

"""
    Kernal functions corresponding to S(ω)
"""
function kernel(n::Integer, τ::Real, ω::Real, β::Real; spectrum::Symbol=:S)
    if spectrum == :S den = 1.0
    elseif spectrum == :A den = 1.0 - exp(-β*ω)
    else @error "spectrum should be :S or :A"
    end
    res =  (-1)^n * exp(-τ*ω) + exp(-(β-τ)*ω)
    res = ω^n * res / 2π
    return res/den
end

function build_kernal(n::Integer, τ::AbstractVector, ω::AbstractVector, β::Real, spectrum::Symbol = :S)
    dω = ω[2] - ω[1]
    K = zeros(length(τ), length(ω))
    for i = 1:length(τ), j = 1:length(ω)
        K[i,j] = kernel(n, τ[i], ω[j], β, spectrum = spectrum)*dω
    end
    return K
end
@non_differentiable build_kernal(n::Integer, τ::AbstractVector, ω::AbstractVector, β::Real)
@non_differentiable build_kernal(n::Integer, τ::AbstractVector, ω::AbstractVector, β::Real, spectrum::Symbol)


"""
    Generate ω range 0:dω:ωmax
"""
function build_range(dx::Number, xmax::Number, xmin::Number=0.)
    res = map(x->x, xmin:dx:xmax)
    len = length(res)
    return res, len
end
@non_differentiable build_range(dx::Number, xmax::Number)


"""
    Fluction-Dissipation Theorem: A(ω) = (1 - e^(-βω)) S(ω)
"""
function ASrelation(ω::Real, S::Real, β::Real)
    return (1.0 - exp(-β*ω)) * S
end

function SArelation(ω::Real, A::Real, β::Real)
    return A / (1.0 - exp(-β*ω))
end


