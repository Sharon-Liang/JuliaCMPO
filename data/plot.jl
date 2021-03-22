using PyPlot
using DelimitedFiles



e = readdlm("./data/0323/fexact-beta-20.txt")
d1 = readdlm("./data/0323/frand-beta-20.txt")
d2 = readdlm("./data/0323/ffix-beta-20.txt")
d3 = readdlm("./data/0323/fhigh-beta-20.txt")

figure # F loss figure
title("TFIM: β=20, χ = 8, date : 2020-03-23")
xlim([0,2])
#ylim([0,1.e-8])
xlabel("Γ/J")
ylabel("F loss")
plot(d1[:,1], d1[:,2] - e[:,2], linewidth = 2)
plot(d2[:,1], d2[:,2] - e[:,2], linewidth = 2, "--")
plot(d3[:,1], d3[:,2] - e[:,2], linewidth = 2, "--")

legend(["random init", "fixed init","high-T init"])

savefig("./figure/Floss-0323-beta-20-1.pdf")
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
d = readdlm("./data/0323/sx-exact-beta-20.txt")
dr = readdlm("./data/0323/sx-rand-beta-20.txt")
df = readdlm("./data/0323/sx-fix-beta-20.txt")
dh = readdlm("./data/0323/sx-high-beta-20.txt")
title("β = 20, χ = 8, date: 2020-03-23")
xlabel("Γ/J")
ylabel("<sx>")
plot(d[:,1], d[:,2], linewidth =2)
plot(dr[:,1], dr[:,2], linewidth =2, "--")
plot(df[:,1], df[:,2], linewidth =2, "--")
plot(dh[:,1], dh[:,2], linewidth =2, "--")
legend(["exact", "random init", "fixed init","high-T init"])

savefig("./figure/sx-0323-beta-20-1.pdf")
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
