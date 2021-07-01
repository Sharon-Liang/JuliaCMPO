using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf

println("2021-06-29: TFIsing-multiple-output.jl")
println("2021-06-29: calculate sz")
χ = 8
x = make_operator(pauli(:x),χ)
z = make_operator(pauli(:z),χ)

gamma = [i for i in range(0.,2.,step = 0.02)]
beta = [1, 20]
mb = maximum(beta)

pcollect1 = Vector{String}(undef, length(beta))
pcollect2 = Vector{String}(undef, length(beta))

for i = 1:length(beta)
    pcollect1[i] = @sprintf "../data_new/b_%i_sz.jld" beta[i]
    pcollect2[i] = @sprintf "../data_new/f_and_sx_b_%i_sz.txt" beta[i]
end 

for path in pcollect1
    jldopen(path, "w") do file
    end
end

for path in pcollect2
    open(path, "w") do file
    end
end

for g in gamma
    key = string(g)
    w = TFIsing(1.0, g, field=:z)
    arr = init_cmps(χ,w) |> toarray
    β = 1
    while β <= mb
        f = arr -> free_energy(arr, w, β)
        gf! = grad_func(f, arr)
        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
        arr = op.minimizer
        res = (minimum(op), arr, Optim.converged(op))
        if res[3] == false println("Not converged Γ = ", key, ", β = ", β) end
        path1 = @sprintf "../data_new/b_%i_sz.jld" β
        path2 = @sprintf "../data_new/f_and_sx_b_%i_sz.txt" β
    
        if path1 in pcollect1
            jldopen(path1, "r+") do file
                write(file, key, res)
            end
        end

        if path2 in pcollect2
            ψ = tocmps(arr)
            #f_exa = free_energy(1.0, g, β)
            #sx_exa = ave_sx(1.0, g, β)
            sz = thermal_average(z,ψ,w,β)
            open(path2, "a") do file
                #writedlm(file,[g f_exa minimum(op) sx_exa sx])
                writedlm(file,[g sz])
            end
        end  
        β += 1
    end
end