"""
using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf
"""

"""Convert time t in seconds to DD:HH:MM:SS"""
function converttime(t::Number)
    (m, s) = divrem(t, 60)
    (h, m) = divrem(m, 60)
    (d, h) = divrem(m, 24)
    if d == 0
        if h == 0
            if m == 0 return @sprintf "time=%.3fs" s
            else return @sprintf "time=%imin %.3f s" m s
            end
        else
            return @sprintf "time=%ih %imin %.3fs" h m s
        end
    else
        return @sprintf "time=%id %ih %imin %.3fs" d h m s
    end
end

println("2021-11-05: TFIsing-change-gamma.jl")
println("χ = 16")

χ = 16
x = make_operator(pauli(:x),χ)

gamma = [0.5, 1.0, 1.5, 2.0]
beta = [i for i in range(1.,40.,step=0.1)]

pcollect1 = Vector{String}(undef, length(gamma))
pcollect2 = Vector{String}(undef, length(gamma))
upcollect = Vector{String}(undef, length(gamma))

for i = 1:length(gamma)
    pcollect1[i] = @sprintf "../data/chi16/g_%.1f.jld" gamma[i]
    pcollect2[i] = @sprintf "../data/chi16/f_and_sx_g_%.1f.txt" gamma[i]
    upcollect[i] = @sprintf "../data/chi16/ug_%.1f.jld" gamma[i]
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

start = time()
for j = 1:length(gamma)
    t1 = time()
    g = gamma[j]
    w = TFIsing(1.0, g)
    arr = init_cmps(χ,w) |> toarray
    path1 = @sprintf "../data/chi16/g_%.1f.jld" g
    path2 = @sprintf "../data/chi16/f_and_sx_g_%.1f.txt" g
    
    #output file path
    upath = @sprintf "../data/chi16/ug_%.1f.jld" g

    for i = 1:length(beta)
        β = beta[i]; key = string(β)
        f = arr -> free_energy(arr, w, β)
        gf! = gradient_function(f)
        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
        arr = op.minimizer
        res = (minimum(op), arr, Optim.converged(op))
        if res[3] == false 
            println("Not converged Γ = ", key, ", β = ", β)
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
        f_exa = free_energy(1.0, g, β)
        f2 = free_energy(ψ2, w, β)
        sx_exa = ave_sx(1.0, g, β)
        sx = thermal_average(x,ψ,w,β)
        open(path2, "a") do file2
            writedlm(file2,[β f_exa minimum(op) f2 sx_exa sx])
        end
    end
    t2 = time()
    println("finish Γ/J = ", g)
    println(converttime(t2 - t1))
end
finish = time()

println("beta in range (", minimum(beta)," , ",maximum(beta),")")
println("Finish!")
println(converttime(finish - start))

