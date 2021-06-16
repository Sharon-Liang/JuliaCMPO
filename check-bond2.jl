using Pkg
Pkg.activate("./")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

println("change gamma for selected beta")
gamma = [i for i in range(0,0.2,length = 100)]
β = 500
path = @sprintf "./data/check-bond2.txt" 
open(path,"w") do file
    for Γ in gamma
	w = TFIsing(1.0, Γ)
        arr = init_cmps(2,w) |> toarray
        f = arr -> free_energy(arr, w, β)
        gf! = grad_func(f, arr)
        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
        arr = op.minimizer
	f_exa = free_energy(1.0, Γ, β)
	loss = minimum(op) - f_exa
        writedlm(file, [Γ loss])
    end
end
