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

# start with high-T states
init_arr = load("./data/0323/sfix-beta-1.jld")

println("χ = ", χ, " ; β = ", β)

open("./data/beta-20-chi-8/S-history.txt","a") do file
    write(file, "\n 2020-3-23 11:30 \n")
    write(file,"unconverged Γ: \n")

    println("Start the high-T init procedure.")
    jldopen("./data/beta-20-chi-8/shigh-beta-20.jld","w") do io
        open("./data/beta-20-chi-8/fhigh-beta-20.txt","a") do io2
            for j in g
                W = TFIsing(1.0, j)
                key = string(j)
                arr = init_arr[key]
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
                write(io, key, arr)
            end
        end
    end
end

println("Finsh calculating states!")
