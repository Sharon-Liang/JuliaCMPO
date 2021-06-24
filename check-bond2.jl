using Pkg
Pkg.activate("./")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

gamma = [i for i in range(0,0.2,length = 100)]
β = 5
path = @sprintf "./data/check-bond2_b%i.txt" β
open(path,"w") do file
    for Γ in gamma
	w = TFIsing(Γ, 1.0)
        arr = init_cmps(2,w) |> toarray
        f = arr -> free_energy(arr, w, β)
        gf! = grad_func(f, arr)
        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
        arr = op.minimizer
	    f_exa = free_energy(Γ, 1.0, β)
	    loss = minimum(op) - f_exa
        writedlm(file, [Γ loss])
    end
end
