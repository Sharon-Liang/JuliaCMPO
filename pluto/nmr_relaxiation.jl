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
beta = [i for i in range(1, 20, step=0.2)]

# ╔═╡ 8bf42a6a-51eb-42b4-b9aa-e6f8a5148d2f
begin
	data = Vector(undef,3)
	for i = 1:3
		data[i] = load(pcollect[i])
	end
	"load data"
end

# ╔═╡ 9dd8334d-7664-44f9-b052-b3858d913066
begin
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
end

# ╔═╡ 05f4d459-5bc7-4dbf-ab13-42041c0cf8e7
begin
	z = make_operator(pauli(:z),8)
	x = make_operator(pauli(:x),8)
	"x, z operator"
end

# ╔═╡ dfe621e5-a860-47e0-b670-2fd4f57eae68
function nmr_relaxiation(g::Real,η::Float64=0.05)
	i = findall(x-> x==g, gamma)[1]
	d = data[i]
	w = TFIsing(1., g)
	len = length(beta)
	res = zeros(len)
	for i = 1:len
		β = beta[i]; key = string(β)
		ψ = tocmps(d[key][2])
		res[i] = structure_factor(0,z,z,ψ,w,β,η=η)
	end
	return res
end

# ╔═╡ 4dac8cf2-8e9d-4aa5-9061-1228332a08b9
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

# ╔═╡ 5f036037-0e4a-42a8-be10-2e467ca1015d
function plot_nmr(g::Real, η::Float64)
	datalab = @sprintf "Γ/J=%.1f,η=%.3f" g η
	t = nmr_relaxiation(g, η)
	if g==1
		xdata = beta
		ydata = t
	else
		xdata = beta[1:6]
		ydata = t[1:6]
	end
	ylimit = (minimum(t)*0.8, maximum(t)*1.1)
	plot(1 ./ beta, t, line=(:dash), marker=(:circle,3),label=datalab)
	plot!(xlabel="T", ylabel = "S(0)")
	f, fname = fitnmr(g, xdata, ydata)
	plot!(1 ./ beta, f(beta), lw = 2, label =fname, ylim = ylimit)
end

# ╔═╡ e72e5cc9-0ad4-42d9-be50-979134a7055a
function plot_nmr(g::Real)
	datalab = @sprintf "Γ/J=%.1f,β/2 average" g
	t = nmr_relaxiation(g)
	xdata = beta
	ydata = t
	ylimit = (minimum(t)*0.8, maximum(t)*1.1)
	plot(1 ./ beta, t, line=(:dash), marker=(:circle,3),label=datalab)
	plot!(xlabel="T", ylabel = "S(0)")
	f, fname = fitnmr(g, xdata, ydata)
	plot!(1 ./ beta, f(beta), lw = 2, label =fname, ylim = ylimit)
end

# ╔═╡ 8c6de546-fc37-411c-8ee8-c8b331056142
plot_nmr(1.0, 0.05)

# ╔═╡ 3645514b-8f48-4607-ad2c-1ed36b15b741
plot_nmr(1.0, 0.1)

# ╔═╡ 69f9a51a-518a-46bc-9cd9-c1be6908f433
plot_nmr(1.0, 0.005)

# ╔═╡ 797747bb-6e76-47a1-8a39-f2ad99158f32
plot_nmr(0.5,0.05)

# ╔═╡ 85b7398c-754d-43ae-93f8-2847c78581e5
plot_nmr(0.5, 0.1)

# ╔═╡ 1bc84f68-08c9-47fa-8e3a-9dc5d7d31d47
plot_nmr(0.5, 0.005)

# ╔═╡ 1729c59d-3c4e-450d-8c68-d742372e0063
plot_nmr(2., 0.05)

# ╔═╡ 734c4847-1f03-478a-9346-90ef3e0a376b
plot_nmr(2.,0.1)

# ╔═╡ 3e95bb33-0c55-45b2-8a16-78a5f1f0fbcc
plot_nmr(2., 0.005)

# ╔═╡ f3b5afce-f8b5-4177-b6cb-86b3e9140972
plot_nmr(2.)

# ╔═╡ 28fb60ad-00c6-4423-8102-fc76c0e8711f
plot_nmr(1.)

# ╔═╡ 812d5668-44bc-4011-8ed5-7e907dc4dcde
plot_nmr(0.5)

# ╔═╡ Cell order:
# ╠═838f933c-db31-4887-b206-562bd060746f
# ╟─f4fda113-17ce-4312-93ff-fb361815b970
# ╟─3f990c4b-97da-4336-8e34-1bd496162328
# ╟─8bf42a6a-51eb-42b4-b9aa-e6f8a5148d2f
# ╟─dfe621e5-a860-47e0-b670-2fd4f57eae68
# ╟─4dac8cf2-8e9d-4aa5-9061-1228332a08b9
# ╟─9dd8334d-7664-44f9-b052-b3858d913066
# ╠═5f036037-0e4a-42a8-be10-2e467ca1015d
# ╠═e72e5cc9-0ad4-42d9-be50-979134a7055a
# ╟─8c6de546-fc37-411c-8ee8-c8b331056142
# ╟─3645514b-8f48-4607-ad2c-1ed36b15b741
# ╟─69f9a51a-518a-46bc-9cd9-c1be6908f433
# ╠═797747bb-6e76-47a1-8a39-f2ad99158f32
# ╟─85b7398c-754d-43ae-93f8-2847c78581e5
# ╟─1bc84f68-08c9-47fa-8e3a-9dc5d7d31d47
# ╟─1729c59d-3c4e-450d-8c68-d742372e0063
# ╟─734c4847-1f03-478a-9346-90ef3e0a376b
# ╟─3e95bb33-0c55-45b2-8a16-78a5f1f0fbcc
# ╠═f3b5afce-f8b5-4177-b6cb-86b3e9140972
# ╠═28fb60ad-00c6-4423-8102-fc76c0e8711f
# ╠═812d5668-44bc-4011-8ed5-7e907dc4dcde
# ╠═05f4d459-5bc7-4dbf-ab13-42041c0cf8e7
# ╠═9f3511fc-e24d-11eb-130d-15fd034ec045
