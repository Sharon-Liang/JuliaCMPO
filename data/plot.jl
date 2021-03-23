using PyPlot
using DelimitedFiles



e1 = readdlm("./data/0323/fexact-gamma-0.1.txt")
e2 = readdlm("./data/0323/fexact-gamma-1.0.txt")
e3 = readdlm("./data/0323/fexact-gamma-2.0.txt")

d1 = readdlm("./data/0323/f-gamma-0.1.txt")
d2 = readdlm("./data/0323/f-gamma-1.0.txt")
d3 = readdlm("./data/0323/f-gamma-2.0.txt")

figure # F loss figure
title("TFIM: χ = 8, date : 2020-03-23")
xlim([1,20])
#ylim([0,1.e-8])
xlabel("β")
ylabel("F loss")
#plot(d1[:,1], d1[:,2] - e1[:,2], linewidth = 2)
#plot(d2[:,1], d2[:,2] - e2[:,2], linewidth = 2)
plot(d3[:,1], d3[:,2] - e3[:,2], linewidth = 2)
ticklabel_format(axis="y", style="scientific",scilimits=(0,0))
legend(["Γ/J = 2.0"])

savefig("./figure/Floss-gamma-2.0.pdf")
PyPlot.display_figs()

figure # F-compare figure
title("TFIM: β=20, χ = 8, date : 2020-03-19")
xlim([0,2])
#ylim([0,1.e-8])
xlabel("Γ/J")
ylabel("F loss: random init - fixed init")
plot(d1[:,1], zeros(length(d1[:,1])),"--k")
plot(d1[:,1], d1[:,2] - d2[:,2], linewidth = 2)

savefig("./figure/Fcompare-0319.pdf")
PyPlot.display_figs()

figure # <s> figure
e1 = readdlm("./data/0323/sx-exact-gamma-0.1.txt")
e2 = readdlm("./data/0323/sx-exact-gamma-1.0.txt")
e3 = readdlm("./data/0323/sx-exact-gamma-2.0.txt")

d1 = readdlm("./data/0323/sx-gamma-0.1.txt")
d2 = readdlm("./data/0323/sx-gamma-1.0.txt")
d3 = readdlm("./data/0323/sx-gamma-2.0.txt")
title("χ = 8, date: 2020-03-23")
xlabel("β")
ylabel("<sx>")

plot(e1[:,1], e1[:,2], linewidth =2)
plot(e2[:,1], e2[:,2], linewidth =2)
plot(e3[:,1], e3[:,2], linewidth =2)
plot(d1[:,1], d1[:,2], linewidth =2, "--")
plot(d2[:,1], d2[:,2], linewidth =2, "--")
plot(d3[:,1], d3[:,2], linewidth =2, "--")
legend(["Γ=0.1 exact", "Γ=1.0 exact", "Γ=2.0 exact",
        "Γ=0.1", "Γ=1.0", "Γ=2.0"])

savefig("./figure/sx-0323.pdf")
PyPlot.display_figs()


figure # <s-s> figure
dr0 = readdlm("./data/beta-20-chi-8/zz-rand-0.txt")
dr1 = readdlm("./data/beta-20-chi-8/zz-rand-1.txt")
dr2 = readdlm("./data/beta-20-chi-8/zz-rand-2.txt")

df0 = readdlm("./data/beta-20-chi-8/zz-fix-0.txt")
df1 = readdlm("./data/beta-20-chi-8/zz-fix-1.txt")
df2 = readdlm("./data/beta-20-chi-8/zz-fix-2.txt")

title("β = 20, χ = 8, date: 2020-03-06")
xlim([0,20])
xlabel("τ")
ylabel("<sz-sz>")
plot(dr0[:,1], dr0[:,2],linewidth = 2)
plot(dr1[:,1], dr1[:,2],linewidth = 2)
plot(dr2[:,1], dr2[:,2],linewidth = 2)

plot(df0[:,1], df0[:,2],linewidth = 2, "--")
plot(df1[:,1], df1[:,2],linewidth = 2, "--")
plot(df2[:,1], df2[:,2],linewidth = 2, "--")
legend(["rand 0.0", "rand 1.0", "rand 2.0", "fix 0.0", "fix 1.0", "fix 2.0"])

savefig("./figure/sz-sz-0306.pdf")
PyPlot.display_figs()
