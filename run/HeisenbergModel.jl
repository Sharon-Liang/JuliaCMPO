using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

χ = 8
Γ = [1.0]
T = [i for i in range(1,0.02,step=0.01)]
println("change T")
for g in Γ
    path = @sprintf "../data/tnew_%.1f.jld" g
    jldopen(path,"w") do file
        w = TFIsing(1.0, g)
        arr = init_cmps(χ,w) |> toarray
        for t in T
            β = 1/t
            key = string(t)
            f = arr -> free_energy(arr, w, β)
            gf! = grad_func(f, arr)
            op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
            arr = op.minimizer
            res = (minimum(op), arr, Optim.converged(op))
            write(file, key, res)
        end
    end
    println("Γ = ", g)
end