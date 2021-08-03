### A Pluto.jl notebook ###
# v0.15.1

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
	path1 = @sprintf "../data/b_%i.jld" β
	path2 = @sprintf "../data/f_and_sx_b_%i.txt" β
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
	path3 = @sprintf "../data/g_%.1f.jld" g
	path4 = @sprintf "../data/f_and_sx_g_%.1f.txt" g
	d3 = load(path3)
	d4 = readdlm(path4)
	string("d3 = ", path3, "; ", "d4 = ", path4)
end

# ╔═╡ 9ac28f13-96cb-4112-9279-e399123ff1b5
begin
	plot(d4[:,1], d4[:,3]-d4[:,2], line=(2), 
		label=@sprintf "Γ/J = %.1f" g)
	plot!(ylabel ="Free energy loss", xlabel = "β",legend=:topleft)
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
	path5 = @sprintf "../data/sz_b_%i.txt" β
	d5 = readdlm(path5)
	string("d5 = ", path5)
end

# ╔═╡ 1f677532-f7b2-463b-b844-586445478e07
begin
	plot(d5[:,1], d5[:,2], line=(:dash), marker=(:circle, 3), 
		label = "β = 20")
	plot!(ylabel ="< σz >", xlabel = "Γ/J",legend=:topright)
end

# ╔═╡ ee2f4a51-b968-4762-9268-ac563c63c644
md"""
##### Entropy and Specific Heat
"""

# ╔═╡ d6cfbcb3-46d0-4dcf-8b12-823533383417
begin
	beta = copy(d4[:,1])
	T = 1 ./ beta
	logT = log.(T)
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
	s_theory = [entropy(1.,1.,b) for b in beta]
	plot( T , zesc_cmpo[3,:], line=(:dash), marker=(:circle, 2,stroke(0)), 
		label = "Γ/J = 1 cmpo")
	plot!(T, s_theory, lw=1, label="Γ/J = 1 theory")
	plot!(ylabel ="entropy", xlabel = "T",legend=:topleft)
end

# ╔═╡ 13f7d046-e5b3-4fd5-8a5e-230e3ee7ad87
begin
	c_theory = [specific_heat(1.,1.,b) for b in beta]
	plot( T , zesc_cmpo[4,:], line=(:dash), marker=(:circle, 3,stroke(0)), 
		label = "Γ/J = 1 cmpo")
	plot!(T, c_theory , lw=2, label="Γ/J = 1 theory")	
	plot!(ylabel ="specific heat", xlabel = "T",legend=:topleft)
end

# ╔═╡ 7f7d5a73-ac56-4b5b-a48d-79ffe3173b82
md"""
##### Partitian Function
"""

# ╔═╡ c6da9c11-8705-4153-8207-67df57d018ea
begin
	plot( T , zesc_cmpo[1,:], line=(:dash), marker=(:circle, 2,stroke(0)), 
		label = "Γ/J = 1 cmpo")
	plot!(ylabel ="partitian", xlabel = "T",legend=:topleft)
end

# ╔═╡ bd424ece-eb25-4e20-a23d-4847ce95c529
z_theory

# ╔═╡ b85869b9-9f4a-40c5-930c-0840b7cfb077
md"""
##### Energy
"""

# ╔═╡ 35e461dd-83f8-404e-b975-03a3d8d6fe46
begin
	e_theory = [energy(1.0,1.0,b) for b in beta]
	plot( T , zesc_cmpo[2,:], line=(:dash), marker=(:circle, 3,stroke(0)), 
		label = "Γ/J = 1 cmpo")
	plot!(T, e_theory, lw=1, label="Γ/J = 1 theory")
	plot!(ylabel ="energy density", xlabel = "T",legend=:topleft)
end

# ╔═╡ 4161b8fb-ad55-475a-a010-5086700f84f2
md"""
##### Spectrum of $\vdash+\dashv$ and $\vdash\dashv$
"""

# ╔═╡ 0c8e0be4-83ee-48b2-b5c7-e22bc597bcb1
num = 10

# ╔═╡ 31c0b5ce-0cbe-4111-ac99-94a0931874b2
begin
	dE_gong = zeros(length(beta), 63); dE_wang = zeros(length(beta), 127)
	for i = 1:length(beta)
		b = beta[i]; key = string(b)
		ψ = tocmps(d3[key][2])
		gong = ψ * ψ
		wang = ψ * w * ψ
		# eigen values and eigen vectors of gong
		e = eigvals(gong)
		dE_gong[i,:] = e[2:end] .- e[1]
		# eigen values and eigen vectors of wang
		e = eigvals(wang)
		dE_wang[i,:] = e[2:end] .- e[1]
	end
	" dE_gong and dE_wang "
end

# ╔═╡ 235e6cc8-2cf2-4f1d-b95d-5219be75e7d2
function shift_level(E::Vector{Float64}; val::Float64=0.5)
	return E .* val/ E[1]
end

# ╔═╡ 3d6f4488-9d52-42e9-aaeb-57eaaa99b3fa
for i = 1:length(beta)
	dE_gong[i,:] = shift_level(dE_gong[i,:])
	dE_wang[i,:] = shift_level(dE_wang[i,:])
end

# ╔═╡ 03085ea0-d3e4-4374-96df-713575e7a082
begin	
	plot(-logT, dE_wang[:,1], 
		line =(:dash, :reds),
		marker=(:circle, 2, :reds, stroke(0)), 
		label = "wang")
	plot!(-logT, dE_wang[:,2:num], 
		line =(:dash, :reds),
		marker=(:circle, :reds, 2, stroke(0)), 
		label = false)
	plot!(-logT, dE_gong[:,1], 
		line =(:dash, :blues),
		marker=(:circle, 2, :blues, stroke(0)), label = "gong")
	plot!(-logT, dE_gong[:,2:10], 
		line =(:dash, :blues),
		marker=(:circle, 2, :blues, stroke(0)), label = false)
	plot!(xlabel="log(β)", yaxis = ((0,5), 0:0.5:10))
end

# ╔═╡ 3fab0dcf-b745-4e7e-b479-12640e4047a5
begin	
	plot(logT, dE_gong[:,1], 
		line =(:dash, :blues),
		marker=(:circle, 2, :blues, stroke(0)), label = "gong")
	plot!(logT, dE_gong[:,2:10], 
		line =(:dash, :blues),
		marker=(:circle, 2, :blues, stroke(0)), label = false)
end

# ╔═╡ Cell order:
# ╟─950b236a-d7da-11eb-1911-4de0c3345bbe
# ╟─9ddc1f98-baf6-4b33-bba5-02fb94c82b8f
# ╟─1b97e39f-dd29-4672-9ea3-457302691e88
# ╟─6ad5cf43-09f6-4988-a9c3-110964fdcbaf
# ╟─583998e3-370a-4d7f-b018-0a1c723067cc
# ╟─94a279c8-3e26-4702-b3a2-beee20d7ef79
# ╟─9ac28f13-96cb-4112-9279-e399123ff1b5
# ╟─09bee69c-5718-433d-9a11-17577f399883
# ╟─f51e9e98-df9e-424b-ab73-22cdd06e6eef
# ╟─b61bce7a-6fbb-419e-bee7-5731c49dcde0
# ╟─1f677532-f7b2-463b-b844-586445478e07
# ╟─ee2f4a51-b968-4762-9268-ac563c63c644
# ╟─d6cfbcb3-46d0-4dcf-8b12-823533383417
# ╟─414e24c4-12a1-4b89-9c37-fd1cacaf0ba5
# ╟─13f7d046-e5b3-4fd5-8a5e-230e3ee7ad87
# ╟─7f7d5a73-ac56-4b5b-a48d-79ffe3173b82
# ╠═c6da9c11-8705-4153-8207-67df57d018ea
# ╠═bd424ece-eb25-4e20-a23d-4847ce95c529
# ╟─b85869b9-9f4a-40c5-930c-0840b7cfb077
# ╟─35e461dd-83f8-404e-b975-03a3d8d6fe46
# ╟─4161b8fb-ad55-475a-a010-5086700f84f2
# ╟─0c8e0be4-83ee-48b2-b5c7-e22bc597bcb1
# ╟─31c0b5ce-0cbe-4111-ac99-94a0931874b2
# ╟─235e6cc8-2cf2-4f1d-b95d-5219be75e7d2
# ╟─3d6f4488-9d52-42e9-aaeb-57eaaa99b3fa
# ╟─03085ea0-d3e4-4374-96df-713575e7a082
# ╟─3fab0dcf-b745-4e7e-b479-12640e4047a5
# ╟─1b076adf-3117-42fd-b1e1-47a9c1c84254
