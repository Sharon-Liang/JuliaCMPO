using cMPO
using Optim

β = 10.
χ = 2

ψ = init_cmps(χ)
arr = toarray(ψ)

W = TFIsing(1.,1.)

FreeEnergy(ψ, W, β)

of(x::Array{Float64, 3}) = OptimFreeEnergy(x::Array{Float64, 3}, W, β)
of!(gx::Array{Float64, 3}, x::Array{Float64,3}) = OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, β)

op = optimize(of, of!, arr, LBFGS())

minimum(op)
