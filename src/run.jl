using cMPO
using Optim
using DelimitedFiles

χ = 10
W = TFIsing(1.,1.)

#len = 100
#β = [i for i in range(0.1,10,length = len)]

len = 50
β = [10^i for i in range(-2,5,length = len)]

open("./data/chi10-log.txt","w") do io
    ψ = init_cmps(χ)
    arr = toarray(ψ)
    for b in β
        of(x::Array{Float64, 3}) = OptimFreeEnergy(x::Array{Float64, 3}, W, b)
        of!(gx::Array{Float64, 3}, x::Array{Float64,3}) = OptimFreeEnergy!(gx::Array{Float64, 3}, x::Array{Float64,3}, W, b)
        op = optimize(of, of!, arr, LBFGS())
        writedlm(io,[b minimum(op)])
        arr = op.minimizer
    end
end
