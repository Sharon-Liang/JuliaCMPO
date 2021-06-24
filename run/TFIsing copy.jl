using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

χ = 8
x = pauli('x')

println("zero J")
beta = [i for i in range(0.1,30,step = 0.5)]

w = TFIsing(0.0, 1.0)
path = @sprintf "../data/j0.jld" 
jldopen(path,"w") do file
    arr = init_cmps(χ, w) |> toarray
    for β in beta
        key = string(β)
        f = arr -> free_energy(arr, w, β)
        gf! = grad_func(f, arr)
        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
        arr = op.minimizer
        res = (minimum(op), arr, Optim.converged(op))
        write(file, key, res)
    end
end

d = load(path)
path2 = @sprintf "../data/f_and_sx_j0.txt" 
open(path2, "w") do file
    for β in beta
        key = string(β); var = β
        ψ = cmps(d[key][2][:,:,1], d[key][2][:,:,2])
        f_exa = free_energy(0., 1., β)
        sx_exa = ave_sx(0., 1., β)
        sx = thermal_average(x,ψ,w,β)
        writedlm(file,[var f_exa d[key][1] sx_exa sx])
    end
end


