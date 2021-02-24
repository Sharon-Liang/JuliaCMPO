using PyPlot
using DelimitedFiles

data = readdlm("./data/chi10.txt")
exact = readdlm("./data/ExactIsing.txt")

title("F - F_exact(J = Gamma = 1)")
xlim([0,10])
plot(data[:,1], data[:,2] - exact[:,2])

savefig("./data/F_difference.pdf")
PyPlot.display_figs()
