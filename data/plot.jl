using PyPlot
using DelimitedFiles

d1 = readdlm("./data/f_and_sx_b_20_rand.txt")
d2 = readdlm("./data/f_and_sx_b_20.txt")
#d3 = readdlm("./data/f_and_sx_g_2.0.txt")

# F loss figure
function plotstyle(phy::String,var::String,var_lim::AbstractVector;
        sty::String="plain")
        title("TFIM : χ = 8, β = 20, date : 2020-03-25")
        xlim(var_lim); xlabel(var)
        ylabel(phy)
        PyPlot.ticklabel_format(axis="y", style=sty,scilimits=(0,0))
end

figure(figsize = [6.,6.])
plot(d1[:,1], d1[:,3] - d1[:,2],
        linewidth = 2, label = "random init")
plot(d1[:,1], d2[:,3] - d2[:,2],
        linewidth = 2, "--", label = "fixed init")
plotstyle("free_energy_loss", "Γ/J", [0,2]; sty = "scientific")
legend()
PyPlot.display_figs()

figure(figsize = [6.,6.])
plot(d1[:,1], d1[:,4], linewidth = 2, label = "exact")
#plot(d1[:,1], d2[:,4], linewidth = 2, label = "β = 20 exact")
#plot(d1[:,1], d3[:,4], linewidth = 2, label = "Γ/J = 2.0 exact")

plot(d1[:,1], d1[:,5], linewidth = 2, label = "random init")
plot(d1[:,1], d2[:,5], linewidth = 2, "--", label = "fixed init")
#plot(d1[:,1], d3[:,5], linewidth = 2, "--", label = "Γ/J = 2.0")
plotstyle("< sx >", "β", [0,2])
legend()
PyPlot.display_figs()
