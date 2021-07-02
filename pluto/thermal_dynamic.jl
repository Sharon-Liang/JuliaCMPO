### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 1b076adf-3117-42fd-b1e1-47a9c1c84254
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

# ╔═╡ 950b236a-d7da-11eb-1911-4de0c3345bbe
md"""
#### Thermal dynamic quantites
##### Free energy loss
"""

# ╔═╡ 9ddc1f98-baf6-4b33-bba5-02fb94c82b8f
β = 20

# ╔═╡ 1b97e39f-dd29-4672-9ea3-457302691e88
begin
	path1 = @sprintf "../data_new/b_%i.jld" β
	path2 = @sprintf "../data_new/f_and_sx_b_%i.txt" β
	d1 = load(path1)
	d2 = readdlm(path2)
	string("d1 = ", path1, "; ", "d2 = ", path2)
end

# ╔═╡ 6ad5cf43-09f6-4988-a9c3-110964fdcbaf
begin
	plot(d2[:,1], d2[:,3]-d2[:,2], line=(:dash), marker=(:circle,3), 
		label=@sprintf "β = %i " β)
	plot!(ylabel ="Free energy loss",  xlabel = "Γ/J",legend=:best)
end

# ╔═╡ 583998e3-370a-4d7f-b018-0a1c723067cc
g = 1.0

# ╔═╡ 94a279c8-3e26-4702-b3a2-beee20d7ef79
begin
	path3 = @sprintf "../data_new/g_%.1f.jld" g
	path4 = @sprintf "../data_new/f_and_sx_g_%.1f.txt" g
	d3 = load(path3)
	d4 = readdlm(path4)
	string("d3 = ", path3, "; ", "d4 = ", path4)
end

# ╔═╡ 9ac28f13-96cb-4112-9279-e399123ff1b5
begin
	plot(d4[:,1], d4[:,3]-d4[:,2], line=(2), 
		label=@sprintf "Γ/J = %.1f" g)
	plot!(ylabel ="Free energy loss", xlabel = "β",legend=:topright)
end

# ╔═╡ 09bee69c-5718-433d-9a11-17577f399883
md"""
##### $\langle \sigma_x \rangle$ and $\langle \sigma_z \rangle$
"""

# ╔═╡ f51e9e98-df9e-424b-ab73-22cdd06e6eef
begin
	plot(d2[:,1], d2[:,5], line=(:dash), marker=(:circle,3), 
		label=@sprintf "β = %i " β)
	plot!(ylabel ="< σx >",  xlabel = "Γ/J",legend=:topleft)
end

# ╔═╡ b61bce7a-6fbb-419e-bee7-5731c49dcde0
begin
	path5 = @sprintf "../data_new/sz_b_%i.txt" β
	d5 = readdlm(path5)
	string("d5 = ", path5)
end

# ╔═╡ 1f677532-f7b2-463b-b844-586445478e07
begin
	plot(d5[:,1], d5[:,2], line=(:dash), marker=(:circle, 3), 
		label = "β = 20")
	plot!(ylabel ="< σz >", xlabel = "Γ/J",legend=:topright)
end

# ╔═╡ d6cfbcb3-46d0-4dcf-8b12-823533383417
begin
	beta = copy(d4[:,1])
	zesc_cmpo = zeros(4, length(beta))
	w = TFIsing(1.0,1.0)
	for i = 1:length(beta)
		b = beta[i]; key = string(b)
		ψ = tocmps(d3[key][2])
		zesc_cmpo[1,i] = partitian!(ψ,w,b).res
		zesc_cmpo[2,i] = energy(ψ,w,b)
		zesc_cmpo[3,i] = entropy(ψ,w,b)
		zesc_cmpo[4,i] = specific_heat(ψ,w,b)
	end
end

# ╔═╡ 414e24c4-12a1-4b89-9c37-fd1cacaf0ba5
begin
	plot(1 ./ beta, zesc_cmpo[4,:], line=(:dash), marker=(:circle, 3), 
		label = "β = 20")
	plot!(ylabel ="Z", xlabel = "Γ/J",legend=:topright)
end

# ╔═╡ 7f7d5a73-ac56-4b5b-a48d-79ffe3173b82
md"""
##### Partitian Function
"""

# ╔═╡ b85869b9-9f4a-40c5-930c-0840b7cfb077
md"""
##### Energy
"""

# ╔═╡ 6a2bf10a-a5e8-468f-8054-5109a2fd9620
md"""
##### Entropy and Specific Heat
"""

# ╔═╡ Cell order:
# ╟─950b236a-d7da-11eb-1911-4de0c3345bbe
# ╟─9ddc1f98-baf6-4b33-bba5-02fb94c82b8f
# ╠═1b97e39f-dd29-4672-9ea3-457302691e88
# ╠═6ad5cf43-09f6-4988-a9c3-110964fdcbaf
# ╟─583998e3-370a-4d7f-b018-0a1c723067cc
# ╟─94a279c8-3e26-4702-b3a2-beee20d7ef79
# ╟─9ac28f13-96cb-4112-9279-e399123ff1b5
# ╟─09bee69c-5718-433d-9a11-17577f399883
# ╟─f51e9e98-df9e-424b-ab73-22cdd06e6eef
# ╟─b61bce7a-6fbb-419e-bee7-5731c49dcde0
# ╟─1f677532-f7b2-463b-b844-586445478e07
# ╠═d6cfbcb3-46d0-4dcf-8b12-823533383417
# ╠═414e24c4-12a1-4b89-9c37-fd1cacaf0ba5
# ╟─7f7d5a73-ac56-4b5b-a48d-79ffe3173b82
# ╟─b85869b9-9f4a-40c5-930c-0840b7cfb077
# ╟─6a2bf10a-a5e8-468f-8054-5109a2fd9620
# ╠═1b076adf-3117-42fd-b1e1-47a9c1c84254
