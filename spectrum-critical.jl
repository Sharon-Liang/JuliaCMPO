### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 8ff529df-40f6-4c3b-9e77-22db44183123
begin
	using Pkg
	Pkg.activate(".")
	using cMPO
	using DelimitedFiles, HDF5, JLD
	using Plots
	using Printf
	using LsqFit
	using LinearAlgebra
	using Arpack
end

# ╔═╡ c9bd8e82-c4e2-11eb-2de5-55b328d8cf9f
md"""
$\Gamma = J = 1$
"""

# ╔═╡ b5d87692-a6e1-4e04-934a-fa1287a93554
#number of eigen states
num = 10

# ╔═╡ 353d724c-1706-4fee-843b-1a128ad27c21
beta = [i for i in range(0.1,30,step = 0.5)]

# ╔═╡ d6acf104-52cd-4c52-8bab-e5e56dd398ca
w = TFIsing(1.0,1.0)

# ╔═╡ fad68feb-09bf-42d9-aa44-660e195e4478
#load data
begin
	data = load("./data/gnew_1.0.jld")
	"data path = ./data/gnew_1.0.jld"
end

# ╔═╡ aa278724-98c5-4e8c-91ee-b890f604699a
function plot_sz(i::Integer, sz::AbstractArray, num::Integer, beta::AbstractVector)
	lab = @sprintf "(%i, %i)" i-1 i-1
	plot(beta, sz[:, i, i], shape=:circle, markersize=3, label = lab)
	for j= i+1 : num
		lab = @sprintf "(%i,%i)" i-1 j-1
		plot!(beta, sz[:,i,j] .* sz[:,j,i], shape=:circle, markersize=3, label = lab)
	end
	tit = @sprintf "Sz(%i,j) ~ β" i-1
	ylab = @sprintf "|Sz(%i,j)|^2" i-1
	plot!(title = tit, xlabel = "β", ylabel = ylab)
end

# ╔═╡ 32aa73f9-0d71-4711-ad07-0c5db5cbf3a4
begin
	z = pauli('z')
	eye = Matrix(1.0I, 8, 8)
	sz = eye ⊗ z ⊗ eye
	"define pauli z"
end

# ╔═╡ 874de82b-5aa6-41c1-8573-1716e6f06fa9
begin
	dE_gong = zeros(length(beta), num-1); dE_wang = zeros(length(beta), num-1)
	dEβ_gong = zeros(length(beta), num-1); dEβ_wang = zeros(length(beta), num-1)
	Sz = zeros(length(beta), num, num)
	for i = 1:length(beta)
		β = beta[i]; key = string(beta[i])
		ψ = tocmps(data[key][2])
		gong = ψ * ψ
		wang = ψ * w * ψ
		# eigen values and eigen vectors of gong
		e, v = eigs(gong, nev=num, which=:SR)
		dE_gong[i,:] = e[2:end] .- e[1]
		dEβ_gong[i,:] = dE_gong[i,:] .* β
		# eigen values and eigen vectors of wang
		e, v = eigs(wang, nev=num, which=:SR)
		dE_wang[i,:] = e[2:end] .- e[1]
		dEβ_wang[i,:] = dE_wang[i,:] .* β
		Sz[i,:,:] = v' * sz * v
	end	
end

# ╔═╡ 570eef1d-36e4-4f9a-a73e-d2c1cc46f1b7
begin 
	title1 = @sprintf "ΔE/T ~ -ln(T)"
	px = log.(beta)
	plot(px, dEβ_gong, ls=:dash, shape=:circle,markersize=2)
	plot!(px, 0.25 .* px .+ 0.7,lw = 1.5, c=:blue,label=("-0.25log(T)"))
	plot!(title= title1, xlabel = "-ln(T)", ylabel = "ΔEm/T",legend=:topleft)
end

# ╔═╡ 2a5c8e15-9db8-4fb7-955b-cce74ba8a00f
begin 
	title2 = @sprintf "ΔE ~ β"
	plot(beta, dE_gong, ls=:dash, shape=:circle,markersize=2)
	plot!(title= title2, xlabel = "β", ylabel = "ΔEm",legend=:topright)
end

# ╔═╡ affce903-876e-43cc-a04b-f141ee95a67b
plot_sz(1,Sz,num,beta)

# ╔═╡ fdd80b7c-9eb6-478d-8756-b40f483c997c
plot_sz(2,Sz,num,beta)

# ╔═╡ 1466a388-55a4-40ed-bebb-797db74922b4
plot_sz(3,Sz,num,beta)

# ╔═╡ 4eb2497b-5ae7-452c-8408-a210d59218ea
plot_sz(4,Sz,num,beta)

# ╔═╡ 177b7b1d-fdbf-4ea5-91c4-3bf506def062
plot_sz(5,Sz,num,beta)

# ╔═╡ 8db5f545-6fc0-4edb-bef0-f3075a1ca669
plot_sz(6,Sz,num,beta)

# ╔═╡ Cell order:
# ╟─c9bd8e82-c4e2-11eb-2de5-55b328d8cf9f
# ╠═b5d87692-a6e1-4e04-934a-fa1287a93554
# ╠═353d724c-1706-4fee-843b-1a128ad27c21
# ╟─d6acf104-52cd-4c52-8bab-e5e56dd398ca
# ╟─fad68feb-09bf-42d9-aa44-660e195e4478
# ╠═874de82b-5aa6-41c1-8573-1716e6f06fa9
# ╠═570eef1d-36e4-4f9a-a73e-d2c1cc46f1b7
# ╠═2a5c8e15-9db8-4fb7-955b-cce74ba8a00f
# ╠═aa278724-98c5-4e8c-91ee-b890f604699a
# ╠═affce903-876e-43cc-a04b-f141ee95a67b
# ╠═fdd80b7c-9eb6-478d-8756-b40f483c997c
# ╠═1466a388-55a4-40ed-bebb-797db74922b4
# ╠═4eb2497b-5ae7-452c-8408-a210d59218ea
# ╠═177b7b1d-fdbf-4ea5-91c4-3bf506def062
# ╠═8db5f545-6fc0-4edb-bef0-f3075a1ca669
# ╠═32aa73f9-0d71-4711-ad07-0c5db5cbf3a4
# ╠═8ff529df-40f6-4c3b-9e77-22db44183123
