using PyPlot
using DelimitedFiles

data = readdlm("./data/chi10-log.txt")
exact = readdlm("./data/ExactIsing-log.txt")

title(" TFIM: J = Gamma = 1, chi = 10")
xlim([10^-2,10^5])
ylim([0.,4*10^-6])
xlabel("beta")
ylabel("F - F_exact")
semilogx(data[:,1], data[:,2] - exact[:,2])

savefig("./data/F_difference-log.pdf")
PyPlot.display_figs()
