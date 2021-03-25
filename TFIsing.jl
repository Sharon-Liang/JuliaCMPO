using cMPO
using Optim

using DelimitedFiles
using JLD, HDF5
using Printf

x = pauli('x')

χ = 8
#gamma = [0.1, 1, 2]
#beta = [i for i in range(1,20,length = 39)]

beta = [10, 20]
gamma = [i for i in range(0,2,length = 41)]

#for Γ in gamma
    #println(" = ", Γ)
    #w = TFIsing(1., Γ)
    #path = @sprintf "./data/g_%.1f.jld" Γ
for β in beta
    println("β = ", β)
    path = @sprintf "./data/b_%i.jld" β
    jldopen(path,"w") do file
        #arr = init_cmps(χ, w) |> toarray
        #for β in beta
        #    key = string(β)
        for Γ in gamma
            key = string(Γ)
            w = TFIsing(1., Γ)
            arr = init_cmps(χ, w) |> toarray
            f = arr -> free_energy(arr, w, β)
            gf! = grad_func(f, arr)
            op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
            arr = op.minimizer
            res = (minimum(op), arr, Optim.converged(op))
            write(file, key, res)
            println(key)
        end
    end

    d = load(path)
    #path2 = @sprintf "./data/f_and_sx_g_%.1f.txt" Γ
    path2 = @sprintf "./data/f_and_sx_b_%i.txt" β
    open(path2, "w") do file
        #for β in beta
        #    key = string(β); val = β
        for Γ in gamma
            w = TFIsing(1., Γ)
            key = string(Γ); var = Γ
            ψ = cmps(d[key][2][:,:,1], d[key][2][:,:,2])
            f_exa = free_energy(1., Γ, β)
            sx_exa = ave_sx(1., Γ, β)
            sx = thermal_average(x,ψ,w,β)
            writedlm(file,[var f_exa d[key][1] sx_exa sx])
        end
    end
end
