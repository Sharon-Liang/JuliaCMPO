using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

println("2021-07-02: TFIsing-dt.jl")
χ = 8
Γ = [1.0]
T = [i for i in range(0.002,0.2,step=0.002)]
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
            if res[3] == false println("Not converged Γ = ", key, ", β = ", β) end
            arr = op.minimizer
            res = (minimum(op), arr, Optim.converged(op))
            write(file, key, res)
        end
    end
    println("Γ = ", g)
end