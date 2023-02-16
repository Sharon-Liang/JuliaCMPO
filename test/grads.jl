using Test, JuliaCMPO
using LinearAlgebra, Parameters
using Random; Random.seed!()
using FiniteDifferences, Zygote

ngradient(f, x) = grad(central_fdm(5, 1), f, x)
zgradient = Zygote.gradient


#=
### *logtrexp* 
=#


