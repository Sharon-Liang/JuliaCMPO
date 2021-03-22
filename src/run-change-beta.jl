using cMPO
using Optim
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5

println("Start calculating free energy!")

χ = 8
Γ = 2.0

len = 40
beta = [i for i in range(1,20,length = len)]

println("Γ = ", Γ)
open("./data/0323/F-history.txt","a") do file
    write(file, "\n 2020-3-23: 16:03 \n")
    write(file,"Free energy calculation unconverged β: \n")
    println("Start the fix init procedure.")
    open("./data/0323/f-gamma-2.0.txt","w") do io
        jldopen("./data/0323/gamma-2.0.jld","w") do io2
            W = TFIsing(1.,Γ)
            ψ = init_cmps(χ, W)
            arr = toarray(ψ)
            for β in beta
                println("β = ", β)
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
                key = string(β)
                arr = op.minimizer
                write(io2, key, arr)
            end
        end
    end
end

println("Finsh calculating free energy!")
