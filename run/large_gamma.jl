using cMPO
using PyPlot
using DelimitedFiles
using JLD, HDF5
using Printf

z = make_operator(pauli(:z),8);
β = 20;
path = @sprintf "./data/b_%i_jchange.jld" β;
d = load(path);


J = 1.; key =string(J);

#apath =@sprintf "./data/a_j_%.1f.txt" J
gpathm =@sprintf "./data/gm_j_%.1f.txt" J
ipath = @sprintf "./data/ig_j_%.1f.txt" J
gpathn =@sprintf "./data/gn_j_%.1f.txt" J
ψ = tocmps(d[key][2]);
w = TFIsing(J,1.);


#G(iωn)
num = 30
freq = [Masubara_freq(i,β,type=:b) for i=1:num]
G = [Masubara_freq_GF(i,z,z,ψ,w,β) for i=1:num]
Gm = [Masubara_freq_T1(i,z,z,ψ,w,β) for i=1:num]
Gn = [1.0im*freq[i]*Masubara_freq_GF(i,z,z,ψ,w,β) for i=1:num]

open(ipath,"w") do file
    for i=1:num
        writedlm(file,[freq[i] real(1.0im*G[i]) imag(1.0im*G[i])])
    end
end

open(gpathm,"w") do file
    for i=1:num
        writedlm(file,[freq[i] real(Gm[i]) imag(Gm[i])])
    end
end

open(gpathn,"w") do file
    for i=1:num
        writedlm(file,[freq[i] real(Gn[i]) imag(Gn[i])])
    end
end

#spectral_density
#omega = [i for i in range(0,5, length=1000)]
#a1 = [spectral_density(J,x,β,η=0.001) for x in omega]
#open(apath,"w") do file
#    for i=1:1000
#        writedlm(file,[omega[i] a1[i]])
#    end
#end

"""
#correlation_2time
tau = [t for t in range(0,β,length=100)]

gt1 = [correlation_2time(t,z,z,ψ,w,β) for t in tau]
gt2 = [correlation_2time(0,t,β) for t in tau]
dg = gt1 .- gt2

gt3 = [correlation_2time(t,z0,z0,ψ0,w,β) for t in tau]
dg2 = gt3 .- gt2

τ = tau ./ β
figure()
plot(τ, gt2, linewidth=2, label="theory")
plot(τ, gt1, "--o", linewidth = 2, label = "χ=8 cmpo")
plot(τ, gt3, "--*", linewidth = 2, label = "χ=2 cmpo")
xlabel("τ/β")
ylabel("G(τ)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_gf.pdf",bbox_inches="tight")
PyPlot.display_figs()

figure()
plot(τ, abs.(dg), linewidth = 2, label = "χ=8 cmpo")
plot(τ, abs.(dg2),linewidth = 2, label = "χ=2 cmpo")
xlabel("τ/β")
ylabel("G(τ) difference (abs)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_gf_diff.pdf",bbox_inches="tight")
PyPlot.display_figs()

#spectral_density
omega = [i for i in range(0,5, length=1000)]
a1 = [spectral_density(1.,x,β,η=0.001) for x in omega]
open("./data/Agamma1.txt","w") do file
    for i=1:1000
        writedlm(file,[omega[i] a1[i]])
    end
end

a2 = [spectral_density(x,z,z, ψ, w,β) for x in omega]
a3 = [spectral_density(x,z0,z0, ψ0, w,β) for x in omega]

figure()
plot(omega, a1, linewidth=2, label="theory")
plot(omega, a2, "--o", linewidth = 2, label = "χ=8 cmpo")
plot(omega, a3, "--*", linewidth = 2, label = "χ=2 cmpo")
xlabel("ω")
ylabel("Imχ(ω)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_imag_chi.pdf",bbox_inches="tight")

figure()
plot(omega, abs.(a2 .- a1), linewidth = 2, label = "χ=8 cmpo")
plot(omega, abs.(a3 .- a1),"--", linewidth = 2, label = "χ=2 cmpo")
xlabel("ω")
ylabel("Imχ(ω) difference (abs)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_imag_chi_diff.pdf",bbox_inches="tight")

#structure_factor
omega = [i for i in range(1.e-5,4, length=100)]
a1 = [structure_factor(0.,x,β) for x in omega]
a2 = [structure_factor(x,z,z, ψ, w,β) for x in omega]
a3 = [structure_factor(x,z0,z0, ψ0, w,β) for x in omega]

figure()
plot(omega, a1, linewidth=2, label="theory")
plot(omega, a2, "--o", linewidth = 2, label = "χ=8 cmpo")
plot(omega, a3, "--*", linewidth = 2, label = "χ=2 cmpo")
xlabel("ω")
ylabel("S(ω)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_struc.pdf",bbox_inches="tight")

figure()
plot(omega, abs.(a2 .- a1), linewidth = 2, label = "χ=8 cmpo")
plot(omega, abs.(a3 .- a1),"--", linewidth = 2, label = "χ=2 cmpo")
xlabel("ω")
ylabel("S(ω) difference (abs)")
legend()
title(@sprintf "β = %i" β)
PyPlot.savefig("atomic_struc_diff.pdf",bbox_inches="tight")
"""