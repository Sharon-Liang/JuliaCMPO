using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

χ = 8
z = pauli('z')

β = 20
path = @sprintf "../data/bnew_%.1f_sz2.jld" β

gamma = [i for i in range(0,3,step = 0.02)]
jldopen(path,"w") do file
    for Γ in gamma
        key = string(Γ)
        w = TFIsing(1.0, Γ; field='z')
        arr = init_cmps(χ,w) |> toarray
        f = arr -> free_energy(arr, w, β)
        gf! = grad_func(f, arr)
        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
        arr = op.minimizer
        res = (minimum(op), arr, Optim.converged(op))
        write(file, key, res)
    end
end

d = load(path)
path2 = @sprintf "../data/sz_bnew_%.1f_s2.txt" β
open(path2, "w") do file
   for Γ in gamma
        key = string(Γ); var = Γ
        ψ = cmps(d[key][2][:,:,1], d[key][2][:,:,2])
        w = TFIsing(1.0, Γ; field='z')
        sz = thermal_average(z,ψ,w,β)
        writedlm(file,[var sz])
    end
end