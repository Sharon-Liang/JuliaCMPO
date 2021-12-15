using Pkg
Pkg.activate("../")
using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using Printf

"""
println("2021-12-09: xxz model ∂ReG_∂ωn")
    D = 8
    N = 40
    #beta = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0]
    beta = [30.0, 40.0]

    model = "xxz"
    folder = "imagtime"
    dir = @sprintf "../data/%s/%s" model folder
    isdir(dir) || mkdir(dir)

    op = [:+, :-, :x, :iy, :z]
    nop= ["pm", "mp", "px", "py", "pz"]

    Delta = [0.0, 1.0]
    for δ = 1:length(Delta)
        Δ = Delta[δ]
        w = XXZmodel(Δ)
        name = @sprintf "Jz_%.1f" Δ
        data_path = @sprintf "../data/xxz/%s_D_%i.jld" name D
        data = load(data_path)
        for o = 1:length(op)
            op_name = nop[o]
            op1 = make_operator(pauli(op[o]), D)
            op2 = make_operator(pauli(op[o]), 2D)

            for b = 1:length(beta)
                β = beta[b]; key = string(β)
                ψ1 = tocmps(data[key][2]); ψ2 = w * ψ1
                    
                func = "dReG"
                dir1 = @sprintf "%s/%s" dir func
                isdir(dir1) || mkdir(dir1)
                path1 = @sprintf "%s/%s_%s_D_%i_beta_%i.txt" dir1 name op_name D β
                path2 = @sprintf "%s/%s_%s_D_%im2_beta_%i.txt" dir1 name op_name D β

                ωn = [Masubara_freq(i,β) for i=1:N]
                Gt1 = [∂ReG_∂ωn(i,op1,op1',ψ1,w,β) for i=1:N]
                Gt2 = [∂ReG_∂ωn(i,op2,op2',ψ2,w,β) for i=1:N]  

                open(path1,"w") do file
                    for i=1: N
                        writedlm(file,[ωn[i] Gt1[i]])
                    end
                end
                
                open(path2,"w") do file
                    for i=1: N
                        writedlm(file,[ωn[i] Gt2[i]])
                    end
                end
            end
        end
        println("finish ", name)
    end
"""

"""
println("2021-12-12: ising model ∂ReG_∂ωn")
    D = 8
    N = 40
    beta = [1.0, 2.0, 4.0, 6.0, 10.0, 20.0, 30.0, 40.0]

    model = "ising"
    folder = "imagtime"
    dir = @sprintf "../data/%s/%s" model folder
    isdir(dir) || mkdir(dir)

    op1 = make_operator(pauli(:z), D)
    op2 = make_operator(pauli(:z), 2D)

    gamma=[1.0]

    for ga = 1:length(gamma)
        g = gamma[ga]
        w = TFIsing(1.0, g)
        data_path = @sprintf "../data/ising/D_%i/g_%.1f.jld" D g
        data = load(data_path)
        for b = 1:length(beta)
            β = beta[b]; key = string(β)
            ψ1 = tocmps(data[key][2]); ψ2 = w * ψ1
                    
            func = "dReG"
            dir1 = @sprintf "%s/%s" dir func
            isdir(dir1) || mkdir(dir1)
            path1 = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir1 g D β
            path2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i.txt" dir1 g D β

            ωn = [Masubara_freq(i,β) for i=1:N]
            Gt1 = [∂ReG_∂ωn(i,op1,op1',ψ1,w,β) for i=1:N]
            Gt2 = [∂ReG_∂ωn(i,op2,op2',ψ2,w,β) for i=1:N]  

            open(path1,"w") do file
                for i=1: N
                    writedlm(file,[ωn[i] Gt1[i]])
                end
            end
                
            open(path2,"w") do file
                for i=1: N
                    writedlm(file,[ωn[i] Gt2[i]])
                end
            end
        end
        println("finish Γ/J = ", g)
    end
"""

println("2021-12-15: ising model S0iwn")
    D = 8
    N = 40
    beta = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0]

    model = "ising"
    folder = "imagtime"
    dir = @sprintf "../data/%s/%s" model folder
    isdir(dir) || mkdir(dir)

    op1 = make_operator(pauli(:z), D)
    op2 = make_operator(pauli(:z), 2D)

    gamma=[1.0]

    for ga = 1:length(gamma)
        g = gamma[ga]
        w = TFIsing(1.0, g)
        data_path = @sprintf "../data/ising/D_%i/g_%.1f.jld" D g
        data = load(data_path)
        for b = 1:length(beta)
            β = beta[b]; key = string(β)
            ψ1 = tocmps(data[key][2]); ψ2 = w * ψ1
                    
            func = "S0iwn"
            dir1 = @sprintf "%s/%s" dir func
            isdir(dir1) || mkdir(dir1)
            path1 = @sprintf "%s/g_%.1f_D_%i_beta_%i.txt" dir1 g D β
            path2 = @sprintf "%s/g_%.1f_D_%im2_beta_%i.txt" dir1 g D β

            ωn = [Masubara_freq(i,β) for i=1:N]
            Gt1 = [structure_factor(0,op1,op1',ψ1,w,β, η=ωn[i]) for i=1:N]
            Gt2 = [structure_factor(0,op2,op2',ψ2,w,β, η=ωn[i]) for i=1:N]  

            open(path1,"w") do file
                for i=1: N
                    writedlm(file,[ωn[i] Gt1[i]])
                end
            end
                
            open(path2,"w") do file
                for i=1: N
                    writedlm(file,[ωn[i] Gt2[i]])
                end
            end
        end
        println("finish Γ/J = ", g)
    end


"""
println("2021-12-09: xxz model imagtime")
    D = 8
    N = 40
    beta = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0]
    len = 1601

    model = "xxz"
    folder = "imagtime"
    dir = @sprintf "../data/%s/%s" model folder
    isdir(dir) || mkdir(dir)


    op = [:+, :-, :x, :iy, :z]
    nop= ["pm", "mp", "px", "py", "pz"]

    omega = [i for i in range(-4π,4π,length=len)]
    eta = [0.05, 0.001]
    Delta=[0.0, 1.0]

    for δ = 1:length(Delta)
        Δ = Delta[δ]
        w = XXZmodel(Δ)
        name = @sprintf "Jz_%.1f" Δ
        data_path = @sprintf "../data/xxz/%s_D_%i.jld" name D
        data = load(data_path)
        for o = 1:length(op)
            op_name = nop[o]
            op1 = make_operator(pauli(op[o]), D)
            op2 = make_operator(pauli(op[o]), 2D)

            for b = 1:length(beta)
                β = beta[b]; key = string(β)
                ψ1 = tocmps(data[key][2]); ψ2 = w * ψ1
                    
                # G(τ)
                func = "gtau"
                dir1 = @sprintf "%s/%s" dir func
                isdir(dir1) || mkdir(dir1)
                path1 = @sprintf "%s/%s_%s_D_%i_beta_%i.txt" dir1 name op_name D β
                path2 = @sprintf "%s/%s_%s_D_%im2_beta_%i.txt" dir1 name op_name D β

                tau = [i for i in range(0, β, length = len)]
                Gt1 = [correlation_2time(i,op1,op1',ψ1,w,β) for i in tau]
                Gt2 = [correlation_2time(i,op2,op2',ψ2,w,β) for i in tau]  

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
                func = "giwn"
                dir1 = @sprintf "%s/%s" dir func
                isdir(dir1) || mkdir(dir1)
                path1 = @sprintf "%s/%s_%s_D_%i_beta_%i.txt" dir1 name op_name D β
                path2 = @sprintf "%s/%s_%s_D_%im2_beta_%i.txt" dir1 name op_name D β

                ωn = [Masubara_freq(i,β) for i=1:N]
                Gt1 = [Masubara_freq_GF(i,op1,op1',ψ1,w,β) for i=1:N]
                Gt2 = [Masubara_freq_GF(i,op2,op2',ψ2,w,β) for i=1:N]  

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
                func = "gdivwn"
                dir1 = @sprintf "%s/%s" dir func
                isdir(dir1) || mkdir(dir1)
                path1 = @sprintf "%s/%s_%s_D_%i_beta_%i.txt" dir1 name op_name D β
                path2 = @sprintf "%s/%s_%s_D_%im2_beta_%i.txt" dir1 name op_name D β

                ωn = [Masubara_freq(i,β) for i=1:N]
                Gt1 = [Masubara_freq_GFdivOmega(i,op1,op1',ψ1,w,β) for i=1:N]
                Gt2 = [Masubara_freq_GFdivOmega(i,op2,op2',ψ2,w,β) for i=1:N]  

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
        println("finish xxz model for ", name)
    end
    println("finish!")
"""

"""
println("2021-12-09: xxz model spectrum")
    D = 8
    N = 40
    beta = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0]
    len = 1601

    model = "xxz"
    folder = "spectrum-direct-ac"
    dir = @sprintf "../data/%s/%s" model folder
    isdir(dir) || mkdir(dir)

    op = [:+, :-, :x, :iy, :z]
    nop= ["pm", "mp", "px", "py", "pz"]

    omega = [i for i in range(-4π,4π,length=len)]
    eta = [0.05, 0.001]
    Delta=[0.0, 1.0]

    for δ = 1:length(Delta)
        Δ = Delta[δ]
        w = XXZmodel(Δ)
        name = @sprintf "Jz_%.1f" Δ
        data_path = @sprintf "../data/xxz/%s_D_%i.jld" name D
        data = load(data_path)
        for o = 1:length(op)
            op_name = nop[o]
            op1 = make_operator(pauli(op[o]), D)
            op2 = make_operator(pauli(op[o]), 2D)

            for b = 1:length(beta), e=1:length(eta)
                η = eta[e]
                β = beta[b]; key = string(β)
                ψ1 = tocmps(data[key][2]); ψ2 = w * ψ1
                    
                # A(ω)
                func = "Aw"
                dir1 = @sprintf "%s/%s" dir func
                isdir(dir1) || mkdir(dir1)
                path1 = @sprintf "%s/%s_%s_D_%i_beta_%i_eta_%.3f.txt" dir1 name op_name D β η
                path2 = @sprintf "%s/%s_%s_D_%im2_beta_%i_eta_%.3f.txt" dir1 name op_name D β η

                Gt1 = [spectral_density(i,op1,op1',ψ1,w,β) for i in omega]
                Gt2 = [spectral_density(i,op2,op2',ψ2,w,β) for i in omega]  

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

                # S(ω)
                func = "Sw"
                dir1 = @sprintf "%s/%s" dir func
                isdir(dir1) || mkdir(dir1)
                path1 = @sprintf "%s/%s_%s_D_%i_beta_%i_eta_%.3f.txt" dir1 name op_name D β η
                path2 = @sprintf "%s/%s_%s_D_%im2_beta_%i_eta_%.3f.txt" dir1 name op_name D β η

                Gt1 = [structure_factor(i,op1,op1',ψ1,w,β) for i in omega]
                Gt2 = [structure_factor(i,op2,op2',ψ2,w,β) for i in omega]  

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
        println("finish xxz model for ", name)
    end
    println("finish!")
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




