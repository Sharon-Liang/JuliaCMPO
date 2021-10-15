"""
using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using Printf
#using PyPlot
"""

println("2021-10-15: scan β and g")

N = 100

z8 = make_operator(pauli(:z), 8)
z16 = make_operator(pauli(:z), 16)
beta = [10.0, 20.0, 30.0, 40.0]
gamma = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]

for j = 1:length(gamma)
    g = gamma[j]; w = TFIsing(1.,g)
    #g = 1.5; w = TFIsing(1.,g)
    dpath = @sprintf "../data/ug_%.1f.jld" g
    d8 = load(dpath)
    for b = 1:length(beta)
        β = beta[b]; key = string(β)
        ψ8 = d8[key][2] |> tocmps; ψ28 = w * ψ8
    
        path1 = @sprintf "../data/gt/gf_t8_g_%.1f_b_%i.txt" g β
        path2 = @sprintf "../data/gt/gf_t28_g_%.1f_b_%i.txt" g β

        tau = [i for i in range(0, β, length = 100)]
        Gt8 = [correlation_2time(i,z8,z8,ψ8,w,β) for i in tau]
        Gt28 = [correlation_2time(i,z16,z16,ψ28,w,β) for i in tau]  

        open(path1,"w") do file
            for i=1:N
                writedlm(file,[tau[i]/β Gt8[i]])
            end
        end
    
        open(path2,"w") do file
            for i=1:N
                writedlm(file,[tau[i]/β Gt28[i]])
            end
        end
    end
    println("finish Γ/J = ", g)
end

    #structure_factor
    #omega = [i for i in range(0,5,length=100)]
    #s8 = [structure_factor(i,z8,z8,ψ8,w,β,method =:S) for i in omega]
    #s28 = [structure_factor(i,z16,z16,ψ28,w,β,method =:S) for i in omega]
    #s16 = [structure_factor(i,z16,z16,ψ16,w,β,method =:S) for i in omega]

    #path3 = @sprintf "./data/Sw_g_%.1f.txt" g
    #open(path3, "w") do file
    #    for i = 1:length(omega)
    #        writedlm(file,[omega[i] s8[i] s28[i] s16[i]])
    #    end
    #end

    #spectral density
    #omega = [i for i in range(-5,5,length=200)]
    #s8 = [spectral_density(i,z8,z8,ψ8,w,β) for i in omega]
    #s28 = [spectral_density(i,z16,z16,ψ28,w,β) for i in omega]
    #s16 = [spectral_density(i,z16,z16,ψ16,w,β) for i in omega]

    #path3 = @sprintf "./data/ImX_g_%.1f.txt" g
    #open(path3, "w") do file
    #    for i = 1:length(omega)
    #        writedlm(file,[omega[i] s8[i] s28[i] s16[i]])
    #    end
    #end

    # Correlation_2time
    #omega = [i for i in range(0,β,length=100)]
    #s8 = [correlation_2time(i,z8,z8,ψ8,w,β) for i in omega]
    #s28 = [correlation_2time(i,z16,z16,ψ28,w,β) for i in omega]
    #s16 = [correlation_2time(i,z16,z16,ψ16,w,β) for i in omega]

    #path3 = @sprintf "./data/Gt_g_%.1f_β_%i.txt" g β
    #open(path3, "w") do file
        #for i = 1:length(omega)
            #writedlm(file,[omega[i] s8[i] s28[i] s16[i]])
        #end
    #end


    """
    # free_energy
    omega = [i for i in range(1,20,length=100)]
    fe = [free_energy(1.,g,i) for i in omega]
    s8 = [free_energy(ψ8,w,i) for i in omega] .- fe
    s28 = [free_energy(ψ28,w,i) for i in omega] .- fe
    s16 = [free_energy(ψ16,w,i) for i in omega] .- fe

    path3 = @sprintf "./data/floss_g_%.1f.txt" g
    open(path3, "w") do file
        for i = 1:length(omega)
            writedlm(file,[omega[i] s8[i] s28[i] s16[i]])
        end
    end
    """

