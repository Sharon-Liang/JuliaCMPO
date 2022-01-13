using Flux
using Flux: train!
using ChainRulesCore

"""
    Kernal functions
"""
function kernel(n::Integer, τ::Real, ω::Real, β::Real)
    res =  (-1)^n * exp(-τ*ω) + exp(-(β-τ)*ω)
    return ω^n * res / 2π
end

function build_kernal(n::Integer, τ::AbstractVector, ω::AbstractVector, β::Real)
    dω = ω[2] - ω[1]
    K = zeros(length(τ), length(ω))
    for i = 1:length(τ), j = 1:length(ω)
        K[i,j] = kernel(n, τ[i], ω[j], β)*dω
    end
    return K
end
@non_differentiable build_kernal(n::Integer, τ::AbstractVector, ω::AbstractVector, β::Real)


"""
    Generate ω range 0:dω:ωmax
"""
function build_range(dx::Number, xmax::Number)
    return [i for i=0:dx:xmax]
end
@non_differentiable build_range(dx::Number, xmax::Number)


"""
    d^nχ(τ)/dτ^n = ∫dω K(n,τ,ω) S(ω)
"""
function dnχ(n::Integer, τ::AbstractVector, β::Real, spectrum::Function; 
                dω::Number=0.01, ωmax::Real=100) 
    ω = build_range(dω, ωmax)
    S = map(spectrum, ω)
    K = build_kernal(n, τ, ω, β)
    return K*S
end


"""
    Neural Network model
"""
function build_NN_model(N::Integer = 5)
    model = Chain(
        Dense(1, N, softplus, initW =Flux.glorot_normal),
        Dense(N, 1, softplus, initW =Flux.glorot_normal))
    parameters = params(model)
    return model, parameters
end


"""
    loss functions
"""
function loss(data::Tuple, β::Real, spectrum::Function)
    len = div(length(data), 2)
    res = 0.0
    for i = 1:len
        ȳ = dnχ(i-1, data[2*i-1], β, spectrum)
        res += Flux.mse(ȳ, data[2*i])
    end
    return res
end

function loss(data::AbstractVector, β::Real, spectrum::Function)
    return map(d->loss(d, β, spectrum), data)
end






