using cMPO
using Optim
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5

println("Start changing gamma procedure!")

χ = 8
β = 20
len = 51
g = [i for i in range(0,2,length = len)]

println("χ = ", χ, " ; β = ", β)

open("./data/beta-20-chi-8/S-history.txt","a") do file
    write(file, "\n 2020-3-23 10:15 \n")
    write(file,"unconverged Γ: \n")

    println("Start the random init procedure.")

    jldopen("./data/beta-20-chi-8/srand-beta-20.jld","w") do io
        open("./data/beta-20-chi-8/frand-beta-20.txt","a") do io2
            ψ = init_cmps(χ)
            arr = toarray(ψ)
            for j in g
                W = TFIsing(1.0, j)
                of(x::Array{Float64, 3}) =
                OptimFreeEnergy(x::Array{Float64, 3}, W, β)
                of!(gx::Array{Float64, 3}, x::Array{Float64,3}) =
                OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, β)
                op = optimize(of, of!, arr, LBFGS(),Optim.Options(iterations = 10000))
                if Optim.converged(op) == false
                    write(file, string(β))
                    write(file, "\n")
                end
                writedlm(io2,[j minimum(op)])
                arr = op.minimizer
                key = string(j)
                write(io, key, arr)
            end
        end
    end

    println("Start the fix init procedure.")
    jldopen("./data/beta-20-chi-8/sfix-beta-20.jld","w") do io
        open("./data/beta-20-chi-8/ffix-beta-20.txt","a") do io2
            for j in g
                W = TFIsing(1.0, j)
                ψ = init_cmps(χ, W)
                arr = toarray(ψ)
                of(x::Array{Float64, 3}) =
                OptimFreeEnergy(x::Array{Float64, 3}, W, β)
                of!(gx::Array{Float64, 3}, x::Array{Float64,3}) =
                OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, β)
                op = optimize(of, of!, arr, LBFGS(),Optim.Options(iterations = 10000))
                if Optim.converged(op) == false
                    write(file, string(β))
                    write(file, "\n")
                end
                writedlm(io2,[j minimum(op)])
                arr = op.minimizer
                key = string(j)
                write(io, key, arr)
            end
        end
    end
end

println("Finsh calculating states!")
