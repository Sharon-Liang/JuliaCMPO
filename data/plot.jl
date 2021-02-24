using PyPlot
using DelimitedFiles

data = readdlm("./data/chi10.txt")
exact = readdlm("./data/ExactIsing.txt")

title(" TFIM: J = Gamma = 1, chi = 10")
xlim([0,10])
xlabel("beta")
ylabel("F - F_exact")
plot(data[:,1], data[:,2] - exact[:,2])

savefig("./data/F_difference.pdf")
PyPlot.display_figs()
