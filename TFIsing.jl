using Pkg
Pkg.activate("./")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

χ = 8
x = pauli('x')
"""
println("change beta for selected gamma")
gamma = [0.1, 0.5, 0.9, 1, 1.1, 2, 4]
beta = [i for i in range(0.1,30,step = 0.5)]
for Γ in gamma
    w = TFIsing(1.0, Γ)
    path = @sprintf "./data/gnew_%.1f.jld" Γ
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
    path2 = @sprintf "./data/f_and_sx_gnew_%.1f.txt" Γ
    open(path2, "w") do file
       for β in beta
            key = string(β); var = β
            ψ = cmps(d[key][2][:,:,1], d[key][2][:,:,2])
            f_exa = free_energy(1., Γ, β)
            sx_exa = ave_sx(1., Γ, β)
            sx = thermal_average(x,ψ,w,β)
            writedlm(file,[var f_exa d[key][1] sx_exa sx])
        end
    end
    println("Γ = ", Γ)
end
"""

#println("change gamma for selected beta")
println("large Γ region")
J = [i for i in range(0.,1.,step = 0.02)]
beta = [0.1, 1, 10, 20]
for β in beta
    path = @sprintf "./data/jnew_%.1f.jld" β
    jldopen(path,"w") do file
        for j in J
            key = string(j)
            w = TFIsing(j, 1.0)
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
    path2 = @sprintf "./data/f_and_sx_jnew_%.1f.txt" β
    open(path2, "w") do file
       for j in J
            key = string(j); var = j
            ψ = cmps(d[key][2][:,:,1], d[key][2][:,:,2])
            w = TFIsing(j, 1.0)
            f_exa = free_energy(j, 1.0, β)
            sx_exa = ave_sx(j, 1.0, β)
            sx = thermal_average(x,ψ,w,β)
            writedlm(file,[var f_exa d[key][1] sx_exa sx])
        end
    end
    println("finish β = ", β)
end
