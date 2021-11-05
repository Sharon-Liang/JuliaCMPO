using Pkg
Pkg.activate("../")
using cMPO
using Optim
using Printf
using Zygote
using GenericLinearAlgebra

# hessian_reverse function is copied from:
# https://github.com/FluxML/Zygote.jl/blob/master/src/lib/grad.jl
"""
    hessian_reverse(f, x)

This should be equivalent to [`hessian(f, x)`](@ref hessian),
but implemented using reverse over reverse mode, all Zygote.
(This is usually much slower, and more likely to find errors.)
"""
hessian_reverse(f, x::AbstractArray) = jacobian(x -> gradient(f, x)[1], x)[1]
hessian_reverse(f, x::Number) = gradient(x -> gradient(f, x)[1], x)[1]


println("2021-11-05:hessian test")
#set cmps bond dimension χ and temperature
χ = 8; β = 1.0

# cmpo of transverse field Ising model
w = TFIsing(1.0, 1.0);

# init cmps
arr = init_cmps(χ,w) |> toarray;

# free_energy function
f1 = x -> free_energy(x, w, β);

# test code begin
hessian(f1, arr)

hessian_reverse(f1, arr)


#gradient works fine
gradient(f1, arr)[1]
