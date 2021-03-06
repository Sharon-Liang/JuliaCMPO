using PyPlot
using DelimitedFiles



e = readdlm("./data/Gamma-1/ExactIsing.txt")
d1 = readdlm("./data/Gamma-1/F-random.txt")
d2 = readdlm("./data/Gamma-1/F-fix.txt")

figure # F loss figure
title("TFIM: J = Γ = 1, χ = 8, date : 2020-03-06")
xlim([0,20])
#ylim([0,1.e-8])
xlabel("β")
ylabel("F loss")
plot(d1[:,1], d1[:,2] - e[:,2], linewidth = 2)
plot(d2[:,1], d2[:,2] - e[:,2], linewidth = 2, "--")

legend(["random init", "fixed init"])

savefig("./figure/Floss-0306.pdf")
PyPlot.display_figs()

figure # F-compare figure
title("TFIM: J = Γ = 1, χ = 8, date : 2020-03-06")
xlim([0,20])
#ylim([0,1.e-8])
xlabel("β")
ylabel("F loss: random init - fixed init")
plot(d1[:,1], zeros(length(d1[:,1])),"--k")
plot(d1[:,1], d1[:,2] - d2[:,2], linewidth = 2)

savefig("./figure/Fcompare-0306.pdf")
PyPlot.display_figs()

figure # <s> figure
dr = readdlm("./data/beta-20-chi-8/sx-rand.txt")
df = readdlm("./data/beta-20-chi-8/sx-fix.txt")
title("β = 20, χ = 8, date: 2020-03-06")
xlabel("Γ/J")
ylabel("<sx>")
plot(dr[:,1], dr[:,2], linewidth =2)
plot(df[:,1], df[:,2], linewidth =2 ,"--")
legend(["random init", "fixed init"])

savefig("./figure/sx-0306.pdf")
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
