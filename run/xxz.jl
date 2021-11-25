using Pkg
Pkg.activate("../")
using cMPO
using Optim
using DelimitedFiles
using JLD, HDF5
using Printf
using TimerOutputs

println("2021-11-25: XY model and AFM Heisenberg model")

chi = [8, 16]
Delta = [0.0, 1.0]

to = TimerOutput()
println("gradient optim")
for i = 1:2
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
                arr = init_cmps(w, χ) |> toarray 
                for i = 1:length(beta)
                    β = beta[i]; key = string(β)
                    f = arr -> free_energy(arr, w, β)
                    gf! = gradient_function(f)
                    op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
                    arr = op.minimizer
                    res = (minimum(op), arr, Optim.converged(op))
                    if Optim.converged(op) == false 
                        println("Not converged ",nstr,"  model, χ=" ,χ, ", β = ", β)
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

to = TimerOutput()
println("hessian optim")
for i = 1:2
    χ = 8
    Δ = Delta[i]
    Δ == 0 ? nstr = "xx" : nstr = "heisenberg"
    w = XXZmodel(Δ)

    path1 = @sprintf "../data/xxz/%s_D_%i.jld" nstr χ
    path2 = @sprintf "../data/xxz/%s_f_D_%i.txt" nstr χ
    upath = @sprintf "../data/xxz/%s_D_%i_u.jld" nstr χ

    @timeit to nstr begin
        v, dv = init_cmps(w, χ) |> tovector
        for i = 1:length(beta)
            β = beta[i]; key = string(β)
            f = v -> free_energy(v,dv, w, β)
            gf! = gradient_function(f)
            hf! = hessian_function(f)
            
            op = optimize(f, gf!, hf!, v, NewtonTrustRegion(), Optim.Options(iterations=10000))
            v = op.minimizer
            res = (minimum(op), v, dv, Optim.converged(op))
            if Optim.converged(op) == false 
                println("Not converged ",nstr,"  model",", β = ", β)
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
            
            ψ = tocmps(v, dv)
            ψ2 = w * ψ
            f2 = free_energy(ψ2, w, β)
            open(path2, "a") do file2
                writedlm(file2,[β minimum(op) f2])
            end
        end
    end
end

println(to)
println("Finish!")