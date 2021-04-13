using cMPO
using Optim

using DelimitedFiles
using JLD, HDF5
using Printf
using TimerOutputs

#x = pauli('x')
const to = TimerOutput()

chi = [8, 16]
beta = [i for i in range(1,20,step = 0.1)]

Γ = 1.0
w = TFIsing(1., Γ)

for χ in chi
    @timeit to string(χ) begin
        path = @sprintf "./data/g_%.1f_X_%i.jld" Γ χ
        jldopen(path,"w") do file
            arr = init_cmps(χ, w) |> toarray
            for β in beta
                key = string(β)
                f = arr -> free_energy(arr, w, β)
                gf! = grad_func(f, arr)
                op = optimize(f,gf!, arr, LBFGS(),Optim.Options(iterations = 10000))
                arr = op.minimizer
                res = (minimum(op), arr, Optim.converged(op))
                write(file, key, res)
                if rem(β,1) == 0 println(key) end
            end
        end
    end
    show(to; allocations = false)

    d = load(path)
    path2 = @sprintf "./data/f_and_sx_g_%.1f_X_%i.txt" Γ χ
    open(path2, "w") do file
        for β in beta
            key = string(β); val = β
            ψ = cmps(d[key][2][:,:,1], d[key][2][:,:,2])
            f_exa = free_energy(1., Γ, β)
            sx_exa = ave_sx(1., Γ, β)
            sx = thermal_average(x,ψ,w,β)
            writedlm(file,[var f_exa d[key][1] sx_exa sx])
        end
    end
end
