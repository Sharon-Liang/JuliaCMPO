### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 2de56e90-ce4b-11eb-38fc-db0a4da44e78
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

# ╔═╡ 5abee561-14f6-4c37-a669-e7bbeafada32
md"""
### Paramagnetic phase
"""

# ╔═╡ d9759ebc-db50-4024-8ece-44e089815064
num = 40

# ╔═╡ bd9f611b-a9df-4475-aad7-d929168f5add
begin
	aly = readdlm("../data/exact/analytic-fixg100_l.txt")
	for i= 3:1001
		aly[:,i] = (aly[:,i] .- aly[:,2])
	end
	"N=100, Analytic data (aly)"
end

# ╔═╡ cd5ec73c-e76c-419e-b72c-958863ef407b
function load_exact(N::Integer, xlab::Char)
	if xlab == 'Γ'
		path = @sprintf "../data/exact/exact-fixj%i.txt" N
		xl = @sprintf "Γ/J"
	elseif xlab == 'J'
		path = @sprintf "../data/exact/exact-fixg%i.txt" N
		xl = @sprintf "J/Γ"
	else 
		@error " xlabel should be 'Γ' or 'J'. "
	end
	d = readdlm(path)
	return d	
end

# ╔═╡ 22a60fd1-9018-40ed-8011-8d45efa32955
begin
	N1 = 10
	ed1 = load_exact(N1,'J')
	for i= 3:129
		ed1[:,i] = (ed1[:,i] .- ed1[:,2]) .* N1
	end
	@sprintf "exact data (ed1): N = %i " N1
end

# ╔═╡ 4caaf3f6-42f8-445d-b129-45bd936a5dba
begin
	N = 16
	ed = load_exact(N,'J')
	for i= 3:129
		ed[:,i] = (ed[:,i] .- ed[:,2]) .* N
	end
	@sprintf "exact data (ed): N = %i " N
end

# ╔═╡ 946b658d-bccc-44a8-be9f-6686a08cc178
begin
	β = 20
	J = [i for i in range(0.,1.,step = 0.02)]
	path = @sprintf "../data/jnew_%.1f.jld" β
	cd = load(path)
	string("cmpo data path (cd)=", path)
end

# ╔═╡ 9ecfc8e7-0cb5-468a-95f1-f1bd8c222f06
begin
	dE_gong = zeros(length(J), 63); dE_wang = zeros(length(J), 127)
	for i = 1:length(J)
		j = J[i]; key = string(j)
		w = TFIsing(j,1.0)
		ψ = tocmps(cd[key][2])
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

# ╔═╡ c8fc0917-f7f2-4f26-b6f1-ed41cb5d80a3
begin
	plot(J, 2 .* abs.(J .- 1.0),lw = 2, lc=:gray, label=("|Δ|"), fmt=:png)
	plot!(J, dE_wang[:,1], 
		line =(:dash, :reds),
		marker=(:circle, 2, :reds, stroke(0)), label = "cmpo",
		)
	plot!(J, dE_wang[:,2:num], 
		line =(:dash, :reds),
		marker=(:circle, 2, :reds, stroke(0)), label = false
		)
	plot!(J, dE_gong[:,1:num], 
		line =(:dash, :red),
		marker=(:circle, 2, :red), label=false,
		)
	#plot!(ed[:,1], ed[:,3:90], 
	#	line =(:blue),
	#	marker=(:hexagon, 2, :blue))
	plot!(ed1[:,1], ed1[:,3], 
		line =(:blues), label = "ED:N=16",
		marker=(:star, 2, :greens))
	plot!(ed1[:,1], ed1[:,4:50], 
		line =(:blues), label = false,
		marker=(:star, 2, :blue))
	plot!(xlabel="J/Γ", ylabel = "Ei- E0", ylim=(0,4), legend=:bottomleft)
end

# ╔═╡ fbd6e688-9d3b-49a7-a966-a7327105d6f0
begin
	plot(J, 2 .* abs.(J .- 1.0),lw = 2, lc=:green, label=("|Δ|"))
	plot!(J, 4 .* abs.(J .- 1.0),lw = 2, lc=:red, label=("2|Δ|"))
	plot!(aly[:,1], aly[:,3], 
		line =(:blues), label = "Analytic:N=100",
		marker=(:star, 2, :greens))
	plot!(aly[:,1], aly[:,4:280], 
		line =(:blues), label = false,
		marker=(:star, 2, :blue))
	plot!(J, dE_wang[:,1], 
		line =(:dash, :reds),
		marker=(:circle, 2, :reds, stroke(0)), label = "cmpo",
		fmt = :png)
	plot!(J, dE_wang[:,2:num], 
		line =(:dash, :reds),
		marker=(:circle, 2, :reds, stroke(0)), label = false
		)
	#plot!(J, dE_gong[:,1:num], 
	#	line =(:dash, :red),
	#	marker=(:circle, 2, :red)
	#	)
	#plot!(ed[:,1], ed[:,3:90], 
	#	line =(:blue),
	#	marker=(:hexagon, 2, :blue))
	plot!(xlabel="J/Γ", ylabel = "Ei- E0", ylim=(0,3), legend=:bottomleft)
end

# ╔═╡ Cell order:
# ╟─5abee561-14f6-4c37-a669-e7bbeafada32
# ╟─d9759ebc-db50-4024-8ece-44e089815064
# ╠═c8fc0917-f7f2-4f26-b6f1-ed41cb5d80a3
# ╠═fbd6e688-9d3b-49a7-a966-a7327105d6f0
# ╟─9ecfc8e7-0cb5-468a-95f1-f1bd8c222f06
# ╟─22a60fd1-9018-40ed-8011-8d45efa32955
# ╟─4caaf3f6-42f8-445d-b129-45bd936a5dba
# ╟─bd9f611b-a9df-4475-aad7-d929168f5add
# ╟─cd5ec73c-e76c-419e-b72c-958863ef407b
# ╟─946b658d-bccc-44a8-be9f-6686a08cc178
# ╟─2de56e90-ce4b-11eb-38fc-db0a4da44e78
