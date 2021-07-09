### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 3c67426c-e062-11eb-173f-456b7a6fa7ad
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

# ╔═╡ ece72ed3-0601-42f2-a729-0b97dc17ce60
gamma = [0.5, 1.0, 2.0]

# ╔═╡ 36d6be9b-52cf-4232-9948-e70a7b759987
begin
	pcollect = Vector{String}(undef, length(gamma))
	for i = 1:length(gamma)
	    pcollect[i] = @sprintf "../data/g_%.1f.jld" gamma[i]
	end
	"generate data path"
end

# ╔═╡ 0f912cdf-d7e8-44e0-b3ae-40293eba2e62
pcollect

# ╔═╡ f3459d59-0dd3-43b6-abf6-d8517c303cfb
beta = [i for i in range(1, 20, step=0.2)]

# ╔═╡ 47179917-820e-45e8-a7bd-7413b39c70c8
begin
	z = make_operator(pauli(:z),8)
	x = make_operator(pauli(:x),8)
	"x, z operator"
end

# ╔═╡ c1f04fac-54c5-47b9-ae34-7474586d0d08
function nmr_relaxiation(g::Real;eta::Float64=0.05)
	i = findall(x-> x==g, gamma)[1]
	d = load(pcollect[i])
	w = TFIsing(1., g)
	len = length(beta)
	res = zeros(len)
	for i = 1:len
		β = beta[i]; key = string(β)
		ψ = tocmps(d[key][2])
		res[i] = structure_factor(0,z,z,ψ,w,β,η=eta)
	end
	return res
end

# ╔═╡ 3c462abf-ab81-4190-b92f-4480c81ead76
function plot_nmr(g::Real,eta::Float64...)
	plot(1 ./ beta, nmr_relaxiation(g, eta = eta[1]), 
		ls=:dash, shape=:circle,markersize=3,
		label = @sprintf "η = %.3f" eta[1])
	for i = 2:length(eta)
		plot!(1 ./ beta, nmr_relaxiation(g, eta = eta[i]), 
			ls=:dash, shape=:circle,markersize=3,
			label = @sprintf "η = %.3f" eta[i])
	end
	plot!(title = @sprintf " Γ/J = %.1f " g)
	plot!(xlabel ="T", ylabel = "1/T1")
end

# ╔═╡ f2dae0c8-3167-4950-9695-959840c69f4c
plot_nmr(0.5, 0.1, 0.05, 0.005)

# ╔═╡ ea3bbe73-9ce1-468e-8fd8-8ae0bbe2a0b4
plot_nmr(1.0, 0.1, 0.05, 0.005)

# ╔═╡ 027a5698-64af-477e-b136-8e1f40633361
plot_nmr(2.0, 0.1, 0.05, 0.005)

# ╔═╡ 62bc87a5-d8ff-4ab0-9aed-d77e58b92792
plot(1 ./ beta, nmr_relaxiation(0.5))

# ╔═╡ 637948e1-c303-43f6-bea3-1279bf03faf7
plot(1 ./ beta, nmr_relaxiation(1.0))

# ╔═╡ 878ec7e7-4ee1-4188-a137-c226d528a7f5
plot(1 ./ beta, nmr_relaxiation(2.0))

# ╔═╡ Cell order:
# ╟─ece72ed3-0601-42f2-a729-0b97dc17ce60
# ╟─36d6be9b-52cf-4232-9948-e70a7b759987
# ╟─0f912cdf-d7e8-44e0-b3ae-40293eba2e62
# ╟─f3459d59-0dd3-43b6-abf6-d8517c303cfb
# ╟─c1f04fac-54c5-47b9-ae34-7474586d0d08
# ╟─3c462abf-ab81-4190-b92f-4480c81ead76
# ╠═62bc87a5-d8ff-4ab0-9aed-d77e58b92792
# ╠═637948e1-c303-43f6-bea3-1279bf03faf7
# ╠═878ec7e7-4ee1-4188-a137-c226d528a7f5
# ╠═f2dae0c8-3167-4950-9695-959840c69f4c
# ╠═ea3bbe73-9ce1-468e-8fd8-8ae0bbe2a0b4
# ╠═027a5698-64af-477e-b136-8e1f40633361
# ╟─47179917-820e-45e8-a7bd-7413b39c70c8
# ╟─3c67426c-e062-11eb-173f-456b7a6fa7ad
