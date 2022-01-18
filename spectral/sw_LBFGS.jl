include("./StructureFactorTrain.jl")
using DelimitedFiles, Printf
using Optim, Zygote
using FluxOptTools
using BSON: @save, @load

g = 1.0
β = 10.0
D = 8
ω = build_range(0.01, 10)

#load χ(τ) and dχ(τ) data
p1 = @sprintf "./data/ising/imagtime/gtau/g_%.1f_D_%i_beta_%i.txt" g D β
p2 = @sprintf "./data/ising/imagtime/dgtau/g_%.1f_D_%i_beta_%i.txt" g D β

d1 = readdlm(p1)
d2 = readdlm(p2)

x1 = d1[:,1]; y1 = d1[:,2]; data1 = [(x1, y1)]
x2 = d2[:,1]; y2 = d2[:,2]; data2 = [(x1, y1, x2, y2)]

path = @sprintf "./spectral/data/try1/g_%.1f_D_%i_beta_%i_init.bson" g D β 
@load path model
s(ω::Real) = model([ω])[1]

pars = params(model)
Zygote.refresh() # currently needed when defining new adjoints
lossfun, gradfun, fg!, p0 = optfuns(()->loss(data1[1], β, s), pars)
res = Optim.optimize(Optim.only_fg!(fg!), p0, LBFGS(),
            Optim.Options(iterations=1000, store_trace=true, show_trace=true, 
            f_tol = 1.e-6, g_tol = 1.e-5, time_limit = 3600))

sw = map(x->s(x), ω)
spc1 = @sprintf "./spectral/data/try1/g_%.1f_D_%i_beta_%i_chi.txt" g D β
dc1 = readdlm(spc1)

using Plots
plot(ω, dc1[:,2], line=(:black,2), label="10^-9")
plot!(ω, sw, line=(:dash,2), label="10^-7")



max_try = 5
for t = 1:max_try
    println("try = ", t)
    println("build a NN model")
    model, pars = build_NN_model(5)
    s(ω::Real) = model([ω])[1]

    dir = @sprintf "./spectral/data/try%i" t
    isdir(dir) || mkdir(dir)

    #save the init model
    mpi = @sprintf "%s/g_%.1f_D_%i_beta_%i_init.bson" dir g D β 
    @save mpi model

    println("chi optim")
    Zygote.refresh() # currently needed when defining new adjoints
    lossfun, gradfun, fg!, p0 = optfuns(()->loss(data1[1], β, s), pars)
    res = Optim.optimize(Optim.only_fg!(fg!), p0, LBFGS(),
                Optim.Options(iterations=3000, store_trace=true, show_trace=true, 
                time_limit = 3600))

    mp = @sprintf "%s/g_%.1f_D_%i_beta_%i_chi.bson" dir g D β 
    @save mp model
    sp = @sprintf "%s/g_%.1f_D_%i_beta_%i_chi.txt" dir g D β
    open(sp, "w") do file
        for i = 1:length(ω)
            writedlm(file, [ω[i] s(ω[i])])
        end
    end

    trace = Optim.f_trace(res)
    ep = @sprintf "%s/ftrace_g_%.1f_D_%i_beta_%i_chi.txt" dir g D β
    open(ep, "w") do file
        writedlm(file, trace)
    end

    println("load init model and use dchi")
    mpi = @sprintf "%s/g_%.1f_D_%i_beta_%i_init.bson" dir g D β 
    @load mpi model
    pars = params(model)
    s(ω::Real) = model([ω])[1]

    Zygote.refresh() # currently needed when defining new adjoints
    lossfun, gradfun, fg!, p0 = optfuns(()->loss(data2[1], β, s), pars)
    res = Optim.optimize(Optim.only_fg!(fg!), p0, LBFGS(),
                Optim.Options(iterations=3000, store_trace=true, show_trace=true,
                time_limit=3600))

    mp = @sprintf "%s/g_%.1f_D_%i_beta_%i_dchi.bson" dir g D β 
    @save mp model
    sp = @sprintf "%s/g_%.1f_D_%i_beta_%i_dchi.txt" dir g D β
    open(sp, "w") do file
        for i = 1:length(ω)
            writedlm(file, [ω[i] s(ω[i])])
        end
    end
                
    trace = Optim.f_trace(res)
    ep = @sprintf "%s/ftrace_g_%.1f_D_%i_beta_%i_dchi.txt" dir g D β
    open(ep, "w") do file
        writedlm(file, trace)
    end


    println("build another NN model and use dchi")
    model, pars = build_NN_model(5)
    s(ω::Real) = model([ω])[1]
    loss(data2[1],β,s)

    #save the init model
    mpi = @sprintf "%s/g_%.1f_D_%i_beta_%i_init_dchi.bson" dir g D β 
    @save mpi model

    Zygote.refresh() # currently needed when defining new adjoints
    @time lossfun, gradfun, fg!, p0 = optfuns(()->loss(data2[1], β, s), pars)
    @time res = Optim.optimize(Optim.only_fg!(fg!), p0, LBFGS(),
                Optim.Options(iterations=3000, store_trace=true, show_trace=true))

    mp = @sprintf "%s/g_%.1f_D_%i_beta_%i.bson" dir g D β 
    @save mp model
    sp = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir g D β
    open(sp, "w") do file
        for i = 1:length(ω)
            writedlm(file, [ω[i] s(ω[i])])
        end
    end
                
    trace = Optim.f_trace(res)
    ep = @sprintf "%s/ftrace_g_%.1f_D_%i_beta_%i.txt" dir g D β
    open(ep, "w") do file
        writedlm(file, trace)
    end
end