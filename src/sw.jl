using Flux
import Flux.softplus
using Flux: train!
using DelimitedFiles, Printf

softplus(x::AbstractArray) = softplus.(x)

"""
    Kernal functions
    function kernel(τ::Real, ω::Real, β::Real)
        # LoadError: DomainError with -1.0:
        #log will only return a complex result if called with a complex argument. Try log(Complex(x)).
        res =  (-1)^n * exp(-τ*ω) + exp(-(β-τ)*ω)
        return ω^n *res / 2π
    end
"""
function K0(τ::Real, ω::Real, β::Real)
    res =  exp(-τ*ω) + exp(-(β-τ)*ω)
    return res / 2π
end

function K1(τ::Real, ω::Real, β::Real)
    res =  - exp(-τ*ω) + exp(-(β-τ)*ω)
    return ω*res/2π
end

"""
    Neural Network model
"""
NN_model = softplus ∘ Dense(5, 1) ∘ Dense(1, 5)
parameters = params(NN_model)
SF(ω::Real) = NN_model([ω])[1]


"""
    d^nχ(τ)/dτ^n = ∫dω K(n,τ,ω) S(ω)
"""
function χ(τ::Real, β::Real; N::Integer=10000, ω_max::Real=100)
    dω = ω_max / N
    ω = [i for i=1:dω:ω_max]
    res = K0.(τ, ω, β) .* SF.(ω) .* dω |> sum
    return res 
end

function dχ(τ::Real, β::Real; N::Integer=10000, ω_max::Real=100)
    dω = ω_max / N
    ω = [i for i=1:dω:ω_max]
    res = K1.(τ, ω, β) .* SF.(ω) .* dω |> sum
    return res 
end

function loss(τ, y)
    ypred = χ.(τ, β)
    Flux.mse(ypred, y)
end

function loss1(τ, y, τ1, dy)
    ypred = χ.(τ, β)
    dypred = dχ.(τ1, β)
    return Flux.mse(ypred, y) + Flux.mse(dypred, dy)
end

g = 1.0
β = 10.0
D = 8
path = @sprintf "./data/ising/imagtime/gtau/g_%.1f_D_%i_beta_%i.txt" g D β
path1 = @sprintf "./data/ising/imagtime/dgtau/g_%.1f_D_%i_beta_%i.txt" g D β

d = readdlm(path)
d1 = readdlm(path1)

x = hcat(d[:,1]...)
y = hcat(d[:,2]...)


#gradient(()->loss(x, y), parameters)


step = 200
err = zeros(step+1); err[1] = loss(x,y)
opt = Descent(0.05)
data = [(x, y)]
for epoch in 1:step
    train!(loss, parameters, data, opt)
    err[epoch+1] = loss(x,y)
    println(epoch, ": ", err[epoch+1] )
end