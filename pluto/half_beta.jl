### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 81c624b4-054b-11ec-052d-2bcf183bf8b0
begin
	using Pkg
	Pkg.activate("..")
	using cMPO
	using DelimitedFiles, HDF5, JLD
	using Plots
	using Printf
	using LsqFit
	using LinearAlgebra
	using Arpack
	"load packages"
end

# ╔═╡ b2c3386d-968e-4084-bed2-9bc31d277e6b
gamma = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.5, 2.0, 4.0, 6.0, 8.0, 10.0]

# ╔═╡ ec752e02-346f-4d8f-a312-4a588e4cb6f0
begin
	pcollect = Vector{String}(undef, length(gamma))
	for i = 1:length(gamma)
	    pcollect[i] = @sprintf "../data_new/g_%.1f_lowT.jld" gamma[i]
	end
	"generate data path"
end

# ╔═╡ a0f956fc-0c79-43df-a7aa-02745656e2ee
begin
	data = Vector(undef,length(gamma))
	for i = 1:length(gamma)
		data[i] = load(pcollect[i])
	end
	"load data"
end

# ╔═╡ 094afe58-31f8-4564-943a-f562ba381e51
T = [i for i in range(0.1, 1.e-4, length = 200)]

# ╔═╡ fe6e2650-cb2b-47b2-b53a-b884bb64604a
beta = 1 ./ T

# ╔═╡ 169969cf-f6be-48c5-969c-b6fc75585d75
function fitnmr(g::Real, xdata::Vector, ydata::Vector)
	Δ = 2.0*(1-g)
	if Δ == 0.
		@. model1(x,p) = p[1] * x^p[2]
		p0 = zeros(2)
		fit1= curve_fit(model1, xdata, ydata, p0)
		fname = @sprintf "%.3f × T^(-%.3f)" fit1.param[1] fit1.param[2]
		f = x -> model1(x, fit1.param)
	else
		@. model(x,p) = exp(p[2]*Δ*x) * p[1] + p[3]
		p0 = zeros(3)
		fit = curve_fit(model, xdata, ydata, p0)
		fname = (@sprintf "%.3f × exp(%.3fΔ/T) + %.3f" fit.param[1] fit.param[2] 				fit.param[3])
		f = x -> model(x, fit.param)
	end
	return f, fname
end

# ╔═╡ 4b7f2d47-3ce7-4b70-bac0-e6fe8c3a38b6
begin
	z = make_operator(pauli(:z),8)
	x = make_operator(pauli(:x),8)
	"x, z operator"
end

# ╔═╡ d3cff197-b0a5-493d-8be2-227d4e6fa0cc
function nmr_relaxiation(g::Real)
	i = findall(x-> x==g, gamma)[1]
	d = data[i]
	w = TFIsing(1., g)
	len = length(beta)
	res = zeros(len)
	for i = 1:len
		β = beta[i]; key = string(β)
		ψ = tocmps(d[key][2])
		res[i] = π/4*β*correlation_2time(β/2,z,z,ψ,w,β)
	end
	return res
end

# ╔═╡ c9eef1cd-0972-4010-a08e-70d86a34f20a
function plot_nmr(g::Real; fit::Bool=false)
	datalab = @sprintf "Γ/J=%.1f,β/2 average" g
	t = nmr_relaxiation(g)
	xdata = beta
	ydata = t
	ylimit = (minimum(t)*0.8, maximum(t)*1.1)
	plot(1 ./ beta, t, line=(:dash), marker=(:circle,3),label=datalab)
	if fit == true
		f, fname = fitnmr(g, xdata, ydata)
		plot!(1 ./ beta, f(beta), lw = 2, label =fname, ylim = ylimit)
	end
	plot!(xlabel="T", ylabel = "S(0)")
end

# ╔═╡ 834312e5-526d-4eb3-bd1b-b63d08f87379
begin
	plot_nmr(10.0, fit = true)
	#add_plot_nmr(6.0)
	#add_plot_nmr(2.0)
	#add_plot_nmr(1.0)
end

# ╔═╡ 54258e93-1cc8-4ad7-bfbf-b413c556400c
function add_plot_nmr(g::Real; fit::Bool=false)
	datalab = @sprintf "Γ/J=%.1f,β/2 average" g
	t = nmr_relaxiation(g)
	xdata = beta
	ydata = t
	ylimit = (minimum(t)*0.8, maximum(t)*1.1)
	#plot!(xlabel="T", ylabel = "S(0)")
	if fit == true
		f, fname = fitnmr(g, xdata, ydata)
		plot!(1 ./ beta, f(beta), lw = 2, label =fname, ylim = ylimit)
	end
	plot!(1 ./ beta, t, line=(:dash), marker=(:circle,3),label=datalab)
end

# ╔═╡ 6cebf258-d598-4b2e-becc-2092e58c0ea2
begin
	plot_nmr(0.1)
	add_plot_nmr(0.5)
	add_plot_nmr(0.9)
	add_plot_nmr(1.0)
end

# ╔═╡ Cell order:
# ╟─b2c3386d-968e-4084-bed2-9bc31d277e6b
# ╟─ec752e02-346f-4d8f-a312-4a588e4cb6f0
# ╟─a0f956fc-0c79-43df-a7aa-02745656e2ee
# ╟─094afe58-31f8-4564-943a-f562ba381e51
# ╟─fe6e2650-cb2b-47b2-b53a-b884bb64604a
# ╟─d3cff197-b0a5-493d-8be2-227d4e6fa0cc
# ╟─169969cf-f6be-48c5-969c-b6fc75585d75
# ╟─c9eef1cd-0972-4010-a08e-70d86a34f20a
# ╟─54258e93-1cc8-4ad7-bfbf-b413c556400c
# ╠═6cebf258-d598-4b2e-becc-2092e58c0ea2
# ╠═834312e5-526d-4eb3-bd1b-b63d08f87379
# ╟─4b7f2d47-3ce7-4b70-bac0-e6fe8c3a38b6
# ╟─81c624b4-054b-11ec-052d-2bcf183bf8b0
