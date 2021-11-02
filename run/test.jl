using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf
using TimerOutputs

"""timing g = 1.0, β = 20"""
to = TimerOutput()

@timeit to "total running time" begin
    χ = 8; β = 20.0
    w = TFIsing(1.0, 1.0)
    arr = init_cmps(χ,w) |> toarray
    f = arr -> free_energy(arr, w, β)
    @timeit to "gradient calculation" g! = gradient_function(f)
    @timeit to "hessian calculation" h! = hessian_function(f)
    @timeit to "LBFGS" op1 = optimize(f, g!, arr, LBFGS(),Optim.Options(iterations = 10000))
    @timeit to "Newton" op2 = optimize(f, g!, h!, arr, NewtonTrustRegion(),Optim.Options(iterations = 10000))
end

