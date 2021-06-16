### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 1464edb4-c9bf-11eb-2370-af329ba2a3e5
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

# ╔═╡ e59c23f7-d0e4-4232-a132-5a6884c4be7b
num = 15

# ╔═╡ e8d842a8-bc62-4d93-9003-5cdf518a25a4
β = 20

# ╔═╡ 48045ce5-9992-4913-96ee-bef48909af1b
J = [i for i in range(0.,1.,step = 0.02)]

# ╔═╡ fb473514-8e99-4a16-a4b6-a3b85fd5279f
begin
	path = @sprintf "./data/jnew_%.1f.jld" β
	d1 = load(path)
	string("data path =", path)
end

# ╔═╡ 282f5af5-a1e7-4a97-b0a7-ec67aa6b34f6
function plot_sz(i::Integer, sz::AbstractArray, num::Integer, beta::AbstractVector)
	lab = @sprintf "(%i, %i)" i-1 i-1
	plot(beta, sz[:, i, i], shape=:circle, markersize=3, label = lab)
	for j= i+1 : num
		lab = @sprintf "(%i,%i)" i-1 j-1
		plot!(beta, sz[:,i,j] .* sz[:,j,i], shape=:circle, markersize=3, label = lab)
	end
	tit = @sprintf "Sz(%i,j) ~ J/Γ" i-1
	ylab = @sprintf "|Sz(%i,j)|^2" i-1
	plot!(title = tit, xlabel = "J/Γ", ylabel = ylab)
end

# ╔═╡ cf50face-261d-49d1-92ef-582cc5eec9cf
begin
	z = pauli('z')
	eye = Matrix(1.0I, 8, 8)
	sz = eye ⊗ z ⊗ eye
	"define pauli z"
end

# ╔═╡ fc5f8272-41d1-429d-a9a7-4544f8647f53
begin
	dE_gong = zeros(length(J), num-1); dE_wang = zeros(length(J), num-1)
	Sz = zeros(length(J), num, num)
	for i = 1:length(J)
		j = J[i]; key = string(j)
		w = TFIsing(j,1.0)
		ψ = tocmps(d1[key][2])
		gong = ψ * ψ
		wang = ψ * w * ψ
		# eigen values and eigen vectors of gong
		e, v = eigs(gong, nev=num, which=:SR)
		dE_gong[i,:] = e[2:end] .- e[1]
		# eigen values and eigen vectors of wang
		e, v = eigs(wang, nev=num, which=:SR)
		dE_wang[i,:] = e[2:end] .- e[1]
		Sz[i,:,:] = v' * sz * v
	end	
end

# ╔═╡ f6d2309d-c76e-4a3a-98d5-4d82bb2cfa22
begin #plot energy level
	title1 = @sprintf "Gap function in paramagnetic phase at β=%.1f" β
	plot(J, dE_wang, ls=:dash, shape=:circle,markersize=2)
	#plot(temp, T3[:,2], ls=:dash, shape=:circle, label="cMPO:η=0.1")
	plot!(J, 2 .* abs.(J .- 1.0),lw = 2, label=("|Δ|"))
	plot!(J, 4 .* abs.(J .- 1.0),lw = 2, label=("2|Δ|"))
	plot!(title= title1, xlabel = "J/Γ", ylabel = "ΔEm",legend=:topright)
end

# ╔═╡ 56fa8ff0-8789-42ad-a057-7045144d4db5
plot_sz(1, Sz, num, J)

# ╔═╡ 5c6bfcf4-2b2e-483b-a4c7-e112f1350cec
plot_sz(2, Sz, num, J)

# ╔═╡ 94408f21-1a38-43b3-8fe1-2b7ccd67ecd5
plot_sz(3, Sz, num, J)

# ╔═╡ d963f236-d590-4843-9dfc-05efb9dad1fd
plot_sz(4, Sz, num, J)

# ╔═╡ 24dd81ce-abb9-4408-a208-a00637bcf7f9
plot_sz(5, Sz, num, J)

# ╔═╡ Cell order:
# ╠═e59c23f7-d0e4-4232-a132-5a6884c4be7b
# ╟─e8d842a8-bc62-4d93-9003-5cdf518a25a4
# ╟─48045ce5-9992-4913-96ee-bef48909af1b
# ╟─fb473514-8e99-4a16-a4b6-a3b85fd5279f
# ╠═fc5f8272-41d1-429d-a9a7-4544f8647f53
# ╟─f6d2309d-c76e-4a3a-98d5-4d82bb2cfa22
# ╠═56fa8ff0-8789-42ad-a057-7045144d4db5
# ╠═5c6bfcf4-2b2e-483b-a4c7-e112f1350cec
# ╠═94408f21-1a38-43b3-8fe1-2b7ccd67ecd5
# ╠═d963f236-d590-4843-9dfc-05efb9dad1fd
# ╠═24dd81ce-abb9-4408-a208-a00637bcf7f9
# ╟─282f5af5-a1e7-4a97-b0a7-ec67aa6b34f6
# ╟─cf50face-261d-49d1-92ef-582cc5eec9cf
# ╟─1464edb4-c9bf-11eb-2370-af329ba2a3e5
