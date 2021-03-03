using PyPlot
using DelimitedFiles


"""
data = readdlm("./data/Gamma-1/F_energy.txt")
exact = readdlm("./data/Gamma-1/ExactIsing.txt")

title("TFIM: J = Gamma = 1, chi = 10")
xlim([0,20])
#ylim([0,1.e-8])
xlabel("beta")
ylabel("F - F_{exact}")
plot(data[:,1], data[:,2] - exact[:,2])

savefig("./data/F_difference.pdf")
PyPlot.display_figs()
"""

#d1 = readdlm("./data/beta-20/Sz-Sz-0.2.txt")
#d2 = readdlm("./data/beta-20/Sz-Sz-1.txt")
#d3 = readdlm("./data/beta-20/Sz-Sz-2.txt")
data = readdlm("./data/beta-1/Sz.txt")
title("beta = 1, chi = 10, applied field: 1.e-5 z")
xlim([0,2])
#ylim([0, 0.26])
xlabel("Gamma/J")
ylabel("Sz")
#plot(d1[:,1], d1[:,2])
#plot(d2[:,1], d2[:,2])
#plot(d3[:,1], d3[:,2])
#legend(["Gamma = 0.2", "Gamma = 1.0", "Gamma = 2.0"])
plot(data[:,1], data[:,2])
savefig("./figure/Sz-beta-1.pdf")
PyPlot.display_figs()
