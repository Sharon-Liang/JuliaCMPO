"""
using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf
"""
to = TimerOutput()

println("2021-11-19: TFIsing-change-gamma.jl")

χ = 32
println("χ = ", χ)
x = make_operator(pauli(:x),χ)

gamma = [1.0, 1.5, 2.0]
beta = [i for i in range(1.,40.,step=0.1)]

pcollect1 = Vector{String}(undef, length(gamma))
pcollect2 = Vector{String}(undef, length(gamma))
upcollect = Vector{String}(undef, length(gamma))

for i = 1:length(gamma)
    pcollect1[i] = @sprintf "../data/chi32/g_%.1f.jld" gamma[i]
    pcollect2[i] = @sprintf "../data/chi32/f_g_%.1f.txt" gamma[i]
    upcollect[i] = @sprintf "../data/chi32/ug_%.1f.jld" gamma[i]
end 

for path in pcollect1
    jldopen(path, "w") do file
    end
end

for path in pcollect2
    open(path, "w") do file
    end
end

for path in upcollect
    jldopen(path, "w") do file
    end
end

for j = 1:length(gamma)
    t1 = time()
    g = gamma[j]
    w = TFIsing(1.0, g)
    arr = init_cmps(χ,w) |> toarray
    path1 = @sprintf "../data/chi32/g_%.1f.jld" g
    path2 = @sprintf "../data/chi32/f_g_%.1f.txt" g
    
    #output file path
    upath = @sprintf "../data/chi32/ug_%.1f.jld" g

    gname = @sprintf "g = %.1f" g
    @timeit to gname begin
        for i = 1:length(beta)
            β = beta[i]; key = string(β)
            f = arr -> free_energy(arr, w, β)
            gf! = gradient_function(f)

            bname = @sprintf "optim β = %.1f" β
            @timeit to bname begin
                op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
            end
            arr = op.minimizer
            res = (minimum(op), arr, Optim.converged(op))
            if res[3] == false 
                println("Not converged Γ = ", g, ", β = ", β)
            end

            jldopen(path1, "r+") do file
                write(file, key, res)
            end

            d, r = divrem(β,10)
            if r == 0
                jldopen(upath, "r+") do file
                    write(file, key, res)
                end
            end
            
            ψ = tocmps(arr)
            ψ2 = w * ψ
            f_exa, err= free_energy(1.0, g, β)
            f2 = free_energy(ψ2, w, β)
            open(path2, "a") do file2
                writedlm(file2,[β f_exa minimum(op) f2 err])
            end
        end
    end
end
println(to)
println("Finish!")


