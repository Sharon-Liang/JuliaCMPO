### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ bb98f3ec-fe48-11eb-0fa0-5d8f4b23b417
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

# ╔═╡ bad9c632-c7d7-4461-b4b6-c489a3b91319
beta = [i for i in range(1,100,step=0.1)]

# ╔═╡ 9e0bcc39-b7c1-417c-94d2-9bd8935f0b08
begin
	data =  load("../data/heisenberg.jld")
	"../data/heisenberg.jld" 
end

# ╔═╡ 7edf6487-cf1e-4ff9-9129-a640f0e43e73
num = 5

# ╔═╡ 76b3b43a-1fe4-4ba5-b814-ee7bc1755cb5
w = HeisenbergModel()

# ╔═╡ 38a5e88f-7fe9-4bb2-9958-ca9406719cde
begin
	dE_gong = zeros(length(beta), num-1)
	dE_wang = zeros(length(beta), num-1)
	#sz = zeros(ComplexF32,length(beta), num, num)
	for i = 1:length(beta)
		j = beta[i]; key = string(j)
		ψ = tocmps(data[key][2])
		gong = ψ * ψ
		wang = ψ * w * ψ
		# eigen values and eigen vectors of gong
		e, v = eigs(gong, nev=num, which=:SR)
		dE_gong[i,:] = real(e[2:end] .- e[1] )
		# eigen values and eigen vectors of wang
		e, v = eigs(wang, nev=num, which=:SR)
		dE_wang[i,:] = real(e[2:end] .- e[1] )
		#sz[i,:,:] = v' * z * v
	end	
end

# ╔═╡ 1991b76a-7207-4217-b273-81cbfd3860a2
begin #plot energy level
	title1 = @sprintf "Gap function"
	plot(beta, dE_wang, ls=:dash, shape=:circle,markersize=2)
	plot!(title= title1, xlabel = "β", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ 1a1509c0-50ef-4e61-8f4b-0b83b32e6e9d
begin
	z = make_operator(pauli(:z),8)
	x = make_operator(pauli(:x),8)
	"x, z operator"
end

# ╔═╡ b1c6f6b4-a907-4cdf-808c-39ba0ad89e3b
function nmr_relaxiation(η::Float64)
	d = data
	w = HeisenbergModel()
	len = length(beta)
	res = zeros(len)
	for i = 1:len
		β = beta[i]; key = string(β)
		ψ = tocmps(d[key][2])
		res[i] = structure_factor(0,z,z,ψ,w,β,η=η)
	end
	return res
end

# ╔═╡ c091b721-c654-4d47-b448-4890762856a3
function nmr_relaxiation()
	d = data
	w = HeisenbergModel()
	len = length(beta)
	res = zeros(len)
	for i = 1:len
		β = beta[i]; key = string(β)
		ψ = tocmps(d[key][2])
		res[i] = π/4*β*correlation_2time(β/2,z,z,ψ,w,β)
	end
	return res
end

# ╔═╡ cb5fbdf3-0a91-4e55-8833-b5d116c1634d
function plot_nmr(η::Float64)
	datalab = @sprintf "Heisenberg model: η=%.3f" η
	t = nmr_relaxiation(η)
	plot(1 ./ beta, t, line=(:dash), marker=(:circle,3),label=datalab)
	plot!(xlabel="T", ylabel = "S(0)")
end

# ╔═╡ 0249eb44-011a-440a-827b-789df5fe33e1
function plot_nmr()
	datalab = @sprintf "Heisenberg model: β/2 average"
	t = nmr_relaxiation()
	plot(1 ./ beta, t, line=(:dash), marker=(:circle,3),label=datalab)
	plot!(xlabel="T", ylabel = "S(0)")
end

# ╔═╡ cfae3287-603f-4638-85f7-3d7ef65b7f9b
plot_nmr(0.05)

# ╔═╡ 51cd4196-2aca-4dee-9b0d-da30ebc76305
plot_nmr(0.1)

# ╔═╡ cab4489c-5ba1-470f-a4b6-efeff6611070
plot_nmr(0.005)

# ╔═╡ 9d5f161f-eae2-4770-8647-6b61bd208b2a
plot_nmr()

# ╔═╡ 58016c22-5c97-4a1f-bff2-173cbaaf00cd
pwd()

# ╔═╡ Cell order:
# ╟─bad9c632-c7d7-4461-b4b6-c489a3b91319
# ╟─9e0bcc39-b7c1-417c-94d2-9bd8935f0b08
# ╠═7edf6487-cf1e-4ff9-9129-a640f0e43e73
# ╟─76b3b43a-1fe4-4ba5-b814-ee7bc1755cb5
# ╠═38a5e88f-7fe9-4bb2-9958-ca9406719cde
# ╠═1991b76a-7207-4217-b273-81cbfd3860a2
# ╠═b1c6f6b4-a907-4cdf-808c-39ba0ad89e3b
# ╠═c091b721-c654-4d47-b448-4890762856a3
# ╠═cb5fbdf3-0a91-4e55-8833-b5d116c1634d
# ╠═0249eb44-011a-440a-827b-789df5fe33e1
# ╠═cfae3287-603f-4638-85f7-3d7ef65b7f9b
# ╠═51cd4196-2aca-4dee-9b0d-da30ebc76305
# ╠═cab4489c-5ba1-470f-a4b6-efeff6611070
# ╠═9d5f161f-eae2-4770-8647-6b61bd208b2a
# ╟─1a1509c0-50ef-4e61-8f4b-0b83b32e6e9d
# ╟─bb98f3ec-fe48-11eb-0fa0-5d8f4b23b417
# ╠═58016c22-5c97-4a1f-bff2-173cbaaf00cd
