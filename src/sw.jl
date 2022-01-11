using Flux
using Flux: train!
using DelimitedFiles, Printf
using ChainRules


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
    Neural Network model
"""
function build_NN_model(N::Integer = 5)
    model = Dense(N, 1, softplus) ∘ Dense(1, N, softplus)
    parameters = params(model)
    return model, parameters
end

model, parameters = build_NN_model()
spectrum_func(ω::Real) = model([ω])[1]


"""
    d^nχ(τ)/dτ^n = ∫dω K(n,τ,ω) S(ω)
"""

@non_differentiable range(a,b)

function χ(τ::AbstractVector, β::Real; N::Integer=10000, ω_max::Real=100) 
    ω = [i for i in range(0,ω_max,length=N)]
    S = map(spectrum_func, ω)
    K = build_kernal(0, τ, ω, β)
    return K*S
end

function loss(τ::AbstractVector, y::AbstractVector, β::Real)
    ypred = χ(τ, β)
    Flux.mse(ypred, y)
end


g = 1.0
β = 10.0
D = 8
path = @sprintf "./data/ising/imagtime/gtau/g_%.1f_D_%i_beta_%i.txt" g D β
path1 = @sprintf "./data/ising/imagtime/dgtau/g_%.1f_D_%i_beta_%i.txt" g D β

d = readdlm(path)
d1 = readdlm(path1)

x = d[:,1]
y = d[:,2]


@time loss(x,y,β)
gradient(()->loss(x, y, β), parameters)



step = 200
err = zeros(step+1); err[1] = loss(x,y)
opt = Descent(0.05)
data = [(x, y)]
for epoch in 1:step
    train!(loss, parameters, data, opt)
    err[epoch+1] = loss(x,y)
    println(epoch, ": ", err[epoch+1] )
end