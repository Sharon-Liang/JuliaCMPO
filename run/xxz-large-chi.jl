using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf
using TimerOutputs

println("2021-11-29: XY model and AFM Heisenberg model")

chi = [16]
Delta = [0.0, 1.0]
beta = [i for i in range(1.,40.,step=0.1)]

to = TimerOutput()
println("gradient optim")
for i = 1:length(Delta)
    Δ = Delta[i]
    Δ == 0 ? nstr = "xx" : nstr = "heisenberg"
    w = XXZmodel(Δ)
    
    len = length(chi)
    pcollect1 = Vector{String}(undef, len)
    pcollect2 = Vector{String}(undef, len)
    upcollect = Vector{String}(undef, len)
    
    for i = 1:len
        pcollect1[i] = @sprintf "../data/xxz/%s_D_%i.jld" nstr chi[i]
        pcollect2[i] = @sprintf "../data/xxz/%s_f_D_%i.txt" nstr chi[i]
        upcollect[i] = @sprintf "../data/xxz/%s_D_%i_u.jld" nstr chi[i]
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

    @timeit to nstr begin
        for d = 1:len
            χ = chi[d]
            path1 = @sprintf "../data/xxz/%s_D_%i.jld" nstr χ
            path2 = @sprintf "../data/xxz/%s_f_D_%i.txt" nstr χ
            upath = @sprintf "../data/xxz/%s_D_%i_u.jld" nstr χ

            chi_name = @sprintf "χ=%i" χ
            @timeit to chi_name begin
                arr = init_cmps(χ, w) |> toarray 
                for i = 1:length(beta)
                    β = beta[i]; key = string(β)
                    f = arr -> free_energy(arr, w, β)
                    gf! = gradient_function(f)
                    op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
                    arr = op.minimizer

                    n = 1
                    while Optim.converged(op) == false && n<=5
                        op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
                        n += 1
                    end

                    arr = op.minimizer
                    res = (minimum(op), arr, Optim.converged(op))
                    if n > 1
                        println(nstr,"  model, χ=" ,χ, ", β = ", β, ", n = ",n)
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
                    f2 = free_energy(ψ2, w, β)
                    open(path2, "a") do file2
                        writedlm(file2,[β minimum(op) f2])
                    end
                end
            end
        end
    end
end
println(to)
println("Finish!")