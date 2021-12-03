using Pkg
Pkg.activate("../")
using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using Printf

println("2021-12-03: ising A(ω) and S(ω)")

D = 8
println("ising D=",D," begin!")
z1 = make_operator(pauli(:z), D)
z2 = make_operator(pauli(:z), 2D)

beta = [10.0, 20.0, 30.0, 40.0]
gamma = [1.0, 1.5, 2.0]
len = 1601
omega = [i for i in range(-4π,4π,length=len)]
eta = [0.05, 0.001]


for j = 1:length(gamma)
    g = gamma[j]; w = TFIsing(1.,g)
    dpath = @sprintf "../data/ising/D_%i/ug_%.1f.jld" D g
    d1 = load(dpath)
    for b = 1:length(beta)
        β = beta[b]; key = string(β)
        ψ1 = d1[key][2] |> tocmps; ψ2 = w * ψ1
        
        for e = 1:length(eta)
            η = eta[e]
            "A(ω)"
            path1 = @sprintf "../data/ising/spectral/Aw_g_%.1f_D_%i_beta_%i_eta_%.3f.txt" g D β η
            path2 = @sprintf "../data/ising/spectral/Aw_g_%.1f_D_2m%i_beta_%i_eta_%.3f.txt" g D β η

            Gt1 = [spectral_density(i,z1,z1,ψ1,w,β,η=η) for i in omega]
            Gt2 = [spectral_density(i,z2,z2,ψ2,w,β, η=η) for i in omega]  

            open(path1,"w") do file
                for i=1:len
                    writedlm(file,[omega[i] Gt1[i]])
                end
            end
    
            open(path2,"w") do file
                for i=1:len
                    writedlm(file,[omega[i] Gt2[i]])
                end
            end

            "S(ω)"
            path1 = @sprintf "../data/ising/spectral/Sw_g_%.1f_D_%i_beta_%i_eta_%.3f.txt" g D β η
            path2 = @sprintf "../data/ising/spectral/Sw_g_%.1f_D_2m%i_beta_%i_eta_%.3f.txt" g D β η

            Gt1 = [structure_factor(i,z1,z1,ψ1,w,β,η=η) for i in omega]
            Gt2 = [structure_factor(i,z2,z2,ψ2,w,β, η=η) for i in omega]  

            open(path1,"w") do file
                for i=1:len
                    writedlm(file,[omega[i] Gt1[i]])
                end
            end
    
            open(path2,"w") do file
                for i=1:len
                    writedlm(file,[omega[i] Gt2[i]])
                end
            end
        end
    end
    println("finish Γ/J = ", g)
end


"""
D = 32
println("ising D=",D," begin!")
z1 = make_operator(pauli(:z), D)
z2 = make_operator(pauli(:z), 2D)

beta = [10.0, 20.0, 30.0, 40.0]
gamma = [1.0, 1.5, 2.0]
for j = 1:length(gamma)
    g = gamma[j]; w = TFIsing(1.,g)
    dpath = @sprintf "../data/ising/D_%i/ug_%.1f.jld" D g
    d1 = load(dpath)
    for b = 1:length(beta)
        β = beta[b]; key = string(β)
        ψ1 = d1[key][2] |> tocmps; ψ2 = w * ψ1
        
        # G(τ)
        path1 = @sprintf "../data/ising/imagtime/gtau_g_%.1f_D_%i_beta_%i.txt" g D β
        path2 = @sprintf "../data/ising/imagtime/gtau_g_%.1f_D_2m%i_beta_%i.txt" g D β

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
        path1 = @sprintf "../data/ising/imagtime/giwn_g_%.1f_D_%i_beta_%i.txt" g D β
        path2 = @sprintf "../data/ising/imagtime/giwn_g_%.1f_D_2m%i_beta_%i.txt" g D β

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
        path1 = @sprintf "../data/ising/imagtime/gdivwn_g_%.1f_D_%i_beta_%i.txt" g D β
        path2 = @sprintf "../data/ising/imagtime/gdivwn_g_%.1f_D_2m%i_beta_%i.txt" g D β

        ωn = [Masubara_freq(i,β) for i=1:N]
        Gt1 = [Masubara_freq_GFdivOmega(i,z1,z1,ψ1,w,β) for i=1:N]
        Gt2 = [Masubara_freq_GFdivOmega(i,z2,z2,ψ2,w,β) for i=1:N]  

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
"""

"""
#==free_energy compare==#
gamma = [1.0, 1.5, 2.0]
for j = 1:length(gamma)
    g = gamma[j]; w = TFIsing(1.,g)
    p1 = @sprintf "./data/floss_g_%.1f.txt" g
    p2 = @sprintf "./data/hessian/f_g_%.1f.txt" g

    d1 = readdlm(p1)
    d2 = readdlm(p2)

    if all(d1[:,1] .== d2[:,1]) 
        d2[:,3] = d2[:,3] .- d2[:,2]
        d2[:,4] = d2[:,4] .- d2[:,2]
        p3 = @sprintf "./data/floss_g_%.1f_new.txt" g
        open(p3, "w") do file
            for i = 1:length(d1[:,1])
                # ω D=8 D=2*8 D=16 D=2*16 D=2*2*8 hD=8 hD=2*8
                writedlm(file,[d1[i,:]' d2[i,3:4]'])
            end
        end
    else
        printf("please check ω")
    end
end
"""





"""
N = 40
beta = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0]
len = 1601

D = 8
println("Heisenberg D=",D," begin!")

op = [:x, :y, :x] #correct gz

dpath = @sprintf "../data/xxz/heisenberg_D_%i.jld" D
d1 = load(dpath)

w = XXZmodel(1.0)
for o = 1:3
    s = 0.5 * pauli(op[o]); sname = string(op[o])
    z1 = make_operator(s, D)
    z2 = make_operator(s, 2D)

    for b = 1:length(beta)
        β = beta[b]; key = string(β)
        ψ1 = tocmps(d1[key][2]); ψ2 = w * ψ1
            
        # G(τ)
        path1 = @sprintf "../data/xxz/imagtime-xxz/gtau_%s_heisenberg_D_%i_beta_%i.txt" sname D β
        path2 = @sprintf "../data/xxz/imagtime-xxz/gtau_%s_heisenberg_D_2m%i_beta_%i.txt" sname D β

        tau = [i for i in range(0, β, length = len)]
        Gt1 = [correlation_2time(i,z1,z1,ψ1,w,β) for i in tau]
        Gt2 = [correlation_2time(i,z2,z2,ψ2,w,β) for i in tau]  

        open(path1,"w") do file
            for i=1:len
                writedlm(file,[tau[i] Gt1[i]])
            end
        end
        
        open(path2,"w") do file
            for i=1:len
                writedlm(file,[tau[i] Gt2[i]])
            end
        end

        # G(ωn)
        path1 = @sprintf "../data/xxz/imagtime-xxz/giwn_%s_heisenberg_D_%i_beta_%i.txt" sname D β
        path2 = @sprintf "../data/xxz/imagtime-xxz/giwn_%s_heisenberg_D_2m%i_beta_%i.txt" sname D β

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
        path1 = @sprintf "../data/xxz/imagtime-xxz/gdivwn_%s_heisenberg_D_%i_beta_%i.txt" sname D β
        path2 = @sprintf "../data/xxz/imagtime-xxz/gdivwn_%s_heisenberg_D_2m%i_beta_%i.txt" sname D β

        ωn = [Masubara_freq(i,β) for i=1:N]
        Gt1 = [Masubara_freq_GFdivOmega(i,z1,z1,ψ1,w,β) for i=1:N]
        Gt2 = [Masubara_freq_GFdivOmega(i,z2,z2,ψ2,w,β) for i=1:N]  

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
end
println("finish!")
"""



