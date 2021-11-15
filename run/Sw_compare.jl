using Pkg
Pkg.activate("../")
using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using Printf

println("2021-11-15: G(τ) and g(iωn)")

D = 16
z1 = make_operator(pauli(:z), D)
z2 = make_operator(pauli(:z), 2D)

N = 40
beta = [10.0, 20.0, 30.0, 40.0]
gamma = [1.0, 1.5, 2.0]

len = 1601

for j = 1:length(gamma)
    g = gamma[j]; w = TFIsing(1.,g)
    dpath = @sprintf "../data/chi16/ug_%.1f.jld" g
    d1 = load(dpath)
    for b = 1:length(beta)
        β = beta[b]; key = string(β)
        ψ1 = d1[key][2] |> tocmps; ψ2 = w * ψ1
        
        # G(τ)
        path1 = @sprintf "../data/imagtime/gtau_g_%.1f_D_%i_beta_%i.txt" g D β
        path2 = @sprintf "../data/imagtime/gtau_g_%.1f_D_2m%i_beta_%i.txt" g D β

        tau = [i for i in range(0, β, length = len)]
        Gt1 = [correlation_2time(i,z1,z1,ψ1,w,β) for i in tau]
        Gt2 = [correlation_2time(i,z2,z2,ψ2,w,β) for i in tau]  

        open(path1,"w") do file
            for i=1:len
                writedlm(file,[tau[i]/β Gt1[i]])
            end
        end
    
        open(path2,"w") do file
            for i=1:len
                writedlm(file,[tau[i]/β Gt2[i]])
            end
        end

        # G(ωn)
        path1 = @sprintf "../data/imagtime/giwn_g_%.1f_D_%i_beta_%i.txt" g D β
        path2 = @sprintf "../data/imagtime/giwn_g_%.1f_D_2m%i_beta_%i.txt" g D β

        ωn = [Masubara_freq(i,β) for i=1:N]
        Gt1 = [Masubara_freq_GF(i,z1,z1,ψ1,w,β) for i=1:N]
        Gt2 = [Masubara_freq_GF(i,z2,z2,ψ2,w,β) for i=1:N]  

        open(path1,"w") do file
            for i=1:N
                writedlm(file,[ωn[i] real(Gt1[i]) imag(Gt1[i])])
            end
        end
    
        open(path2,"w") do file
            for i=1:N
                writedlm(file,[ωn[i] real(Gt2[i]) imag(Gt2[i])])
            end
        end

        # χ(ωn)/ωn
        path1 = @sprintf "../data/imagtime/gdivwn_g_%.1f_D_%i_beta_%i.txt" g D β
        path2 = @sprintf "../data/imagtime/gdivwn_g_%.1f_D_2m%i_beta_%i.txt" g D β

        ωn = [Masubara_freq(i,β) for i=1:N]
        Gt1 = [Masubara_freq_T1(i,z1,z1,ψ1,w,β) for i=1:N]
        Gt2 = [Masubara_freq_T1(i,z2,z2,ψ2,w,β) for i=1:N]  

        open(path1,"w") do file
            for i=1:N
                writedlm(file,[ωn[i] real(Gt1[i]) imag(Gt1[i])])
            end
        end
    
        open(path2,"w") do file
            for i=1:N
                writedlm(file,[ωn[i] real(Gt2[i]) imag(Gt2[i])])
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
    #free_energy compare
    p1 = @sprintf "../data/f_and_sx_g_%.1f.txt" g
    p2 = @sprintf "../data/chi16/f_and_sx_g_%.1f.txt" g

    d1 = readdlm(p1)
    d2 = readdlm(p2)

    if all(d1[:,1] .== d2[:,1]) 
        for   
    b = [i for i in range(1,40,length=200)]
    f0 = [free_energy(1.,g,i) for i in omega]
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

