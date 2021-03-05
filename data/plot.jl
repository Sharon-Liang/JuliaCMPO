using PyPlot
using DelimitedFiles


"""
data = readdlm("./data/Gamma-1/F_energy1.txt")
exact = readdlm("./data/Gamma-1/ExactIsing.txt")

title("TFIM: J = Γ = 1, χ = 8")
xlim([0,20])
#ylim([0,1.e-8])
xlabel("β")
ylabel("F - F_{exact}")
plot(data[:,1], data[:,2] - exact[:,2])

savefig("./data/F_difference1.pdf")
PyPlot.display_figs()
"""


#d1 = readdlm("./data/beta-20-chi-8/Sz-Sz-0.2.txt")
#d2 = readdlm("./data/beta-20-chi-8/Sz-Sz-1.txt")
#d3 = readdlm("./data/beta-20-chi-8/Sz-Sz-2.txt")
#data1 = readdlm("./data/beta-20/Sx.txt")
data2 = readdlm("./data/beta-20-chi-8/Sx-1.txt")

title("β = 20, χ = 8")
xlim([0,2])
#ylim([0, 0.26])
#xlabel("Γ/J")
ylabel("Sx")
#xlabel("τ")
#ylabel("Sz-Sz")
#plot(d1[:,1], d1[:,2])
#plot(d2[:,1], d2[:,2])
#plot(d3[:,1], d3[:,2])
#legend(["Gamma = 0.2", "Gamma = 1.0", "Gamma = 2.0"])
#plot(data1[:,1], data1[:,2])
plot(data2[:,1],data2[:,2])
#legend(["χ = 10", "χ=8"])
savefig("./figure/newinit-8/Sx-1.pdf")
PyPlot.display_figs()
