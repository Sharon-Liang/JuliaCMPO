using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

println("2021-07-06: HeisenbergModel.jl")
χ = 8
beta = [i for i in range(1,100,step=0.1)]
path = @sprintf "../data/heisenberg.jld"
jldopen(path,"w") do file
    w = HeisenbergModel()
    arr = init_cmps(χ,w) |> toarray
    for β in beta
        key = string(β)
        f = arr -> free_energy(arr, w, β)
        gf! = gradient_function(f, arr)
        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
        arr = op.minimizer
        res = (minimum(op), arr, Optim.converged(op))
        if res[3] == false println("Not converged β = ", β) end
        write(file, key, res)
        if β%10 == 0 println("Finish β = ", β) end
    end
end
println("Finish!")