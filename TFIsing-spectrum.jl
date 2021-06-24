using cMPO
using DelimitedFiles, HDF5, JLD
using PyPlot
using Printf
using LsqFit

z = pauli('z')

d1 = load("./data/bnew_20.0_r.jld")
#d2 = load("./data/g_0.1.jld")
#d3 = load("./data/g_2.0.jld")

J = 1.0; Γ = 1.0; W = TFIsing(J,Γ)
beta = [i for i in range(1,20,length = 39)]
t0 = [NMR_relaxation(J,Γ,x) for x in beta]
t1, t2, t3 = zeros(39), zeros(39), zeros(39)
for i = 1:39
    β = beta[i]; key = string(beta[i])
    ψ = cmps(d1[key][2][:,:,1], d1[key][2][:,:,2])
    t1[i] = NMR_relaxation(z,z,ψ,W,β,η = 0.001)
    t2[i] = NMR_relaxation(z,z,ψ,W,β,η = 0.01)
    t3[i] = NMR_relaxation(z,z,ψ,W,β,η = 0.1)
end

#data fit
@. model1(t,p) = p[1] * t^(3/4)
p01 = [2.13]
@. model(t,p) = p[1] * t^(p[2])
p0 = [2.13,3/4]

fit1 = curve_fit(model, beta[3:end], t1[3:end], p0)
p1 = fit1.param
fit11 = curve_fit(model1, beta[3:end], t1[3:end], p01)
p11 = fit11.param[1]

fit2 = curve_fit(model, beta, t2, p0)
p2 = fit2.param
fit2 = curve_fit(model1, beta, t2, p01)
p21 = fit2.param[1]

figure()
plot(1 ./ beta, t2, "o",markersize = 3, label = "numerical")
plot(1 ./ beta, p21 .* beta .^(3/4) , label = "~T^(-3/4)")
plot(1 ./ beta, p2[1] .* beta .^(p2[2]) , label = "~T^(-1.34)")

xlabel("T")
ylabel("NMR 1/T1")
title("Γ = J = 1.0, η = 0.01, date : 2020-04-13")
legend()
PyPlot.display_figs()

eta = [10^x for x in range(-7,-1,length = 100)]
β1 = 5.0; key1 = string(β1)
ψ1 = cmps(d1[key1][2][:,:,1], d1[key1][2][:,:,2])
t = [NMR_relaxation(z,z,ψ1,W,β1,η = x) for x in eta]
figure()
plot(eta, NMR_relaxation(J,Γ,β1)*ones(100),"k--")
plot(eta,t,"-o",markersize = 3)
xscale("log")
xlabel("η")
ylabel("NMR 1/T1")
title("Γ = J = 1.0 , β = 5, date : 2020-04-13")
PyPlot.display_figs()



"""
J = 1.0; Γ = 1.0; W = TFIsing(J,Γ)
β = 20.0 ; key = string(β)
ψ = cmps(d1[key][2][:,:,1], d1[key][2][:,:,2])

w = [10^x for x in range(-5,1,length = 1000)]
c = [critical_zz_sus(x,β) for x in w]
c0 = [susceptibility(x,z,z,ψ,W,β,η = 0.05) for x in w]
c1 = [imag_susceptibility(x,z,z,ψ,W,β,η = 0.1) for x in w]

figure()
plot(w, c0, "o",markersize = 3, label = "η = 0.05")
plot(w, c1, "o",markersize = 3, label = "η = 0.1")
plot(, c, "k", label = "theoretical")
xscale("log")
xlabel("ω")
ylabel("Im χ(ω)")
title("β = 20.0, Γ = J = 1.0 , date : 2020-04-13")
legend()
PyPlot.display_figs()
"""
