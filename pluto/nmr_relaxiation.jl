### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 9f3511fc-e24d-11eb-130d-15fd034ec045
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

# ╔═╡ 838f933c-db31-4887-b206-562bd060746f
#gamma = [0.1, 0.5, 1.0, 2.0, 4.0, 8.0]
gamma = [0.5, 1.0, 2.0]

# ╔═╡ f4fda113-17ce-4312-93ff-fb361815b970
begin
	pcollect = Vector{String}(undef, length(gamma))
	for i = 1:length(gamma)
	    pcollect[i] = @sprintf "../data/g_%.1f.jld" gamma[i]
	end
	"generate data path"
end

# ╔═╡ 3f990c4b-97da-4336-8e34-1bd496162328
lowT = [i for i in range(10, 20, step=0.04)]

# ╔═╡ c30979b4-97dc-4252-ad32-eb8ead672c8a
highT = [i for i in range(1, 20, step=0.2)] 

# ╔═╡ 8bf42a6a-51eb-42b4-b9aa-e6f8a5148d2f
begin
	data = Vector(undef,length(gamma))
	for i = 1:length(gamma)
		data[i] = load(pcollect[i])
	end
	"load data"
end

# ╔═╡ 9dd8334d-7664-44f9-b052-b3858d913066
function linearfit(xdata::Vector, ydata::Vector)
	@. model(x,p) = p[1] * x + p[2] 
	p0 = zeros(2)
	fit= curve_fit(model, xdata, ydata, p0)
	fname = @sprintf "%.3f x + %.3f" fit.param[1] fit.param[2]
	f = x -> model(x, fit.param)
	return f, fname
end

# ╔═╡ 05f4d459-5bc7-4dbf-ab13-42041c0cf8e7
begin
	z = make_operator(pauli(:z),8)
	x = make_operator(pauli(:x),8)
	"x, z operator"
end

# ╔═╡ dfe621e5-a860-47e0-b670-2fd4f57eae68
function nmr_relaxiation(g::Real,beta::Vector,η::Float64)
	i = findall(x-> x==g, gamma)[1]
	d = data[i]
	w = TFIsing(1., g)
	len = length(beta)
	res = zeros(len)
	for i = 1:len
		β = beta[i]; key = string(β)
		ψ = tocmps(d[key][2])
		res[i] = structure_factor(0,z,z,ψ,w,β,η=η,method=:S)
	end
	return res
end

# ╔═╡ 4dac8cf2-8e9d-4aa5-9061-1228332a08b9
function nmr_relaxiation(g::Real,beta::Vector)
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

# ╔═╡ 5f036037-0e4a-42a8-be10-2e467ca1015d
function plot_nmr(g::Real,beta::Vector, η::Float64; fit::Bool=false)
	datalab = @sprintf "Γ/J=%.1f,η=%.3f" g η
	t = nmr_relaxiation(g, beta, η)
	ydata = log.(t)
	if g==1
		xdata = log.(1 ./beta)
		ydata = -ydata
		xlab = "logT"
		ylab = "-log(S0)"

	else
		xdata = beta
		xlab = "β"
		ylab = "log(S0)"
	end
	plot(xdata, ydata, line=(:dash), marker=(:circle,3),label=datalab)
		
	if fit == true
		f, fname = linearfit(xdata, ydata)
		plot!(xdata, f(xdata), lw = 2, label=fname)
	end
	plot!(xlabel=xlab, ylabel = ylab)
end

# ╔═╡ 5dfd4e3b-59ef-499a-b8cf-78f008c3a92c
function add_plot_nmr(g::Real,beta::Vector, η::Float64; fit::Bool=false)
	datalab = @sprintf "Γ/J=%.1f,η=%.3f" g η
	t = nmr_relaxiation(g,beta, η)
	ydata = log.(t)
	if g==1
		xdata = log.(1 ./beta)
		ydata = -ydata
		xlab = "logT"
		ylab = "-log(S0)"

	else
		xdata = beta
		xlab = "β"
		ylab = "log(S0)"
	end
	plot!(xdata, ydata, line=(:dash), marker=(:circle,3),label=datalab)
		
	if fit == true
		f, fname = linearfit(xdata, ydata)
		plot!(xdata, f(xdata), lw = 2, label=fname)
	end
	plot!(xlabel=xlab, ylabel = ylab)
end

# ╔═╡ e72e5cc9-0ad4-42d9-be50-979134a7055a
function plot_nmr(g::Real,beta::Vector; fit::Bool=false)
	datalab = @sprintf "Γ/J=%.1f,β/2 average" g
	t = nmr_relaxiation(g,beta)
	ydata = log.(t)
	if g==1
		xdata = log.(1 ./beta)
		ydata = -ydata
		xlab = "logT"
		ylab = "-log(S0)"

	else
		xdata = beta
		xlab = "β"
		ylab = "log(S0)"
	end
	plot(xdata, ydata, line=(:dash), marker=(:circle,3),label=datalab)
		
	if fit == true
		f, fname = linearfit(xdata, ydata)
		plot!(xdata, f(xdata), lw = 2, label=fname)
	end
	plot!(xlabel=xlab, ylabel = ylab)
end

# ╔═╡ 22b7795f-5109-4c1d-acaf-dcb128d3022e
function add_plot_nmr(g::Real,beta::Vector;fit::Bool=false)
	datalab = @sprintf "Γ/J=%.1f,β/2 average" g
	t = nmr_relaxiation(g, beta)
	ydata = log.(t)
	if g==1
		xdata = log.(1 ./beta)
		ydata = -ydata
		xlab = "logT"
		ylab = "-log(S0)"
	else
		xdata = beta
		xlab = "β"
		ylab = "log(S0)"
	end
	plot!(xdata, ydata, line=(:dash), marker=(:circle,3),label=datalab)
		
	if fit == true
		f, fname = linearfit(xdata, ydata)
		plot!(xdata, f(xdata), lw = 2, label=fname)
	end
	plot!(xlabel=xlab, ylabel = ylab)
end

# ╔═╡ fe18cbe8-bc63-4af0-bba1-52f9242373b4
begin
	plot_nmr(1.0, lowT, 0.05, fit=true)
	add_plot_nmr(1.0, lowT, 0.1, fit=true)
	add_plot_nmr(1.0, lowT, fit = true)
	plot!(legend=:topleft)
end

# ╔═╡ e6968b92-f0f7-4667-8dc4-dff6941450ed
begin
	plot_nmr(0.5, lowT, 0.05, fit=true)
	add_plot_nmr(0.5, lowT, 0.1, fit=true)
	add_plot_nmr(0.5, lowT, fit = true)
	plot!(legend=:topleft)
end

# ╔═╡ 049e4f9a-1312-4a5e-b4a6-b1f2a21c452d
begin
	plot_nmr(0.5, highT, 0.05)
	add_plot_nmr(0.5, highT, 0.1)
	add_plot_nmr(0.5, highT)
	plot!(legend=:bottomright)
end

# ╔═╡ be7908c2-1faf-4f2b-b369-a803842af8a5
begin
	plot_nmr(2.0, lowT, 0.05, fit=true)
	add_plot_nmr(2.0, lowT, 0.1, fit=true)
	add_plot_nmr(2.0, lowT, fit = true)
	plot!(legend=:bottomleft)
end

# ╔═╡ 7e90ccea-8268-4839-bd95-bdc818567d65
begin
	plot_nmr(2.0, highT, 0.05)
	add_plot_nmr(2.0, highT, 0.1)
	add_plot_nmr(2.0, highT)
	plot!(legend=:bottomleft)
end

# ╔═╡ 6b8708f3-6504-423b-b962-ceb69f0551a3
pwd()

# ╔═╡ Cell order:
# ╠═838f933c-db31-4887-b206-562bd060746f
# ╟─f4fda113-17ce-4312-93ff-fb361815b970
# ╟─3f990c4b-97da-4336-8e34-1bd496162328
# ╟─c30979b4-97dc-4252-ad32-eb8ead672c8a
# ╠═8bf42a6a-51eb-42b4-b9aa-e6f8a5148d2f
# ╟─dfe621e5-a860-47e0-b670-2fd4f57eae68
# ╟─4dac8cf2-8e9d-4aa5-9061-1228332a08b9
# ╟─9dd8334d-7664-44f9-b052-b3858d913066
# ╟─5f036037-0e4a-42a8-be10-2e467ca1015d
# ╟─5dfd4e3b-59ef-499a-b8cf-78f008c3a92c
# ╟─e72e5cc9-0ad4-42d9-be50-979134a7055a
# ╟─22b7795f-5109-4c1d-acaf-dcb128d3022e
# ╠═fe18cbe8-bc63-4af0-bba1-52f9242373b4
# ╠═e6968b92-f0f7-4667-8dc4-dff6941450ed
# ╠═049e4f9a-1312-4a5e-b4a6-b1f2a21c452d
# ╠═be7908c2-1faf-4f2b-b369-a803842af8a5
# ╠═7e90ccea-8268-4839-bd95-bdc818567d65
# ╟─05f4d459-5bc7-4dbf-ab13-42041c0cf8e7
# ╟─9f3511fc-e24d-11eb-130d-15fd034ec045
# ╟─6b8708f3-6504-423b-b962-ceb69f0551a3
