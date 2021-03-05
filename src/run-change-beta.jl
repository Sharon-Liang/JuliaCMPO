using cMPO
using Optim
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5

println("Start calculating free energy!")

χ = 8

len = 100
beta = [i for i in range(0.2,20,length = len)]

open("./data/Gamma-1/F-history.txt","a") do file
    write(file, "\n 2020-3-5 \n")
    write(file,"Free energy calculation unconverged β: \n")

    println("Start the random init procedure.")

    write(file, "random init with χ = 8, β = 20 \n")
    open("./data/Gamma-1/F-random.txt","w") do io
        W = TFIsing(1.,1.)
        ψ = init_cmps(χ)
        arr = toarray(ψ)
        for β in beta
            of(x::Array{Float64, 3}) =
            OptimFreeEnergy(x::Array{Float64, 3}, W, β)
            of!(gx::Array{Float64, 3}, x::Array{Float64,3}) =
            OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, β)
            op = optimize(of, of!, arr, LBFGS(),Optim.Options(iterations = 10000))
            if Optim.converged(op) == false
                write(file, string(β))
                write(file, "\n")
            end
            writedlm(io,[β minimum(op)])
            arr = op.minimizer
        end
    end

    println("Start the fix init procedure.")
    write(file, "fix init with χ = 8, β = 20 \n")
    open("./data/Gamma-1/F-fix.txt","w") do io
        W = TFIsing(1.,1.)
        ψ = init_cmps(χ, W)
        arr = toarray(ψ)
        for β in beta
            of(x::Array{Float64, 3}) =
            OptimFreeEnergy(x::Array{Float64, 3}, W, β)
            of!(gx::Array{Float64, 3}, x::Array{Float64,3}) =
            OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, β)
            op = optimize(of, of!, arr, LBFGS(),Optim.Options(iterations = 10000))
            if Optim.converged(op) == false
                write(file, string(β))
                write(file, "\n")
            end
            writedlm(io,[β minimum(op)])
            arr = op.minimizer
        end
    end
end

println("Finsh calculating free energy!")
