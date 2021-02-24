using cMPO
using Optim
using DelimitedFiles

χ = 10
ψ = init_cmps(χ)
arr = toarray(ψ)

W = TFIsing(1.,1.)

len = 10
β = [i for i in range(20,30,length = len)]

open("./data/lowT.txt","w") do io
    for b in β
        of(x::Array{Float64, 3}) = OptimFreeEnergy(x::Array{Float64, 3}, W, b)
        of!(gx::Array{Float64, 3}, x::Array{Float64,3}) = OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, b)
        op = optimize(of, of!, arr, LBFGS())
        writedlm(io,[b minimum(op)])
    end
end
