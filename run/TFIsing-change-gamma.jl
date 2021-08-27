"""
using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf
"""

println("2021-08-25: TFIsing-change-gamma.jl")

χ = 8
x = make_operator(pauli(:x),χ)

gamma = [0.1, 0.3, 0.7, 4.0, 6.0, 8.0, 10.0]
beta = [i for i in range(10,1000,step=0.02)]
#T = [i for i in range(0.1, 1.e-4, length = 200)]
#beta = 1 ./ T
mg = maximum(gamma)

pcollect1 = Vector{String}(undef, length(gamma))
pcollect2 = Vector{String}(undef, length(gamma))

for i = 1:length(gamma)
    pcollect1[i] = @sprintf "../data/g_%.1f.jld" gamma[i]
    pcollect2[i] = @sprintf "../data/f_and_sx_g_%.1f.txt" gamma[i]
end 

for path in pcollect1
    jldopen(path, "w") do file
    end
end

for path in pcollect2
    open(path, "w") do file
    end
end

for j = 1:length(gamma)
    g = gamma[j]
    w = TFIsing(1.0, g)
    arr = init_cmps(χ,w) |> toarray
    path1 = @sprintf "../data/g_%.1f.jld" g
    path2 = @sprintf "../data/f_and_sx_g_%.1f.txt" g

    for i = 1:length(beta)
        β = beta[i]; key = string(β)
        f = arr -> free_energy(arr, w, β)
        gf! = gradient_function(f, arr)
        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
        arr = op.minimizer
        res = (minimum(op), arr, Optim.converged(op))
        if res[3] == false 
            println("Not converged Γ = ", key, ", β = ", β)
        end

        jldopen(path1, "r+") do file
            write(file, key, res)
        end
        
        ψ = tocmps(arr)
        f_exa = free_energy(1.0, g, β)
        sx_exa = ave_sx(1.0, g, β)
        sx = thermal_average(x,ψ,w,β)
        open(path2, "a") do file2
            writedlm(file2,[β f_exa minimum(op) sx_exa sx])
        end
    end
    println("finish Γ/J = ", g)
end

println("Finish!")