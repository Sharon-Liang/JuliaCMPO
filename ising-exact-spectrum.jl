### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 66817958-c4da-11eb-04ec-5f000da4c04c
begin
	using DelimitedFiles
	using JLD, HDF5
	using Plots
end

# ╔═╡ 27c49980-570a-4845-9d7a-7c540f0d2abf
begin
	g7 = readdlm("./data/exact/exact-fixg7.txt")
	for i= 3:129
		g7[:,i] = g7[:,i] .- g7[:,2]
	end
	plot(g7[:,1], g7[:,3:end],ls=:dash, shape=:circle,markersize=2)
	plot!(title= "N=7", xlabel = "J/Γ", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ 166f0878-1ed4-4ca1-8118-fb650d65acce
begin
	g8 = readdlm("./data/exact/exact-fixg8.txt")
	for i= 3:129
		g8[:,i] = g8[:,i] .- g8[:,2]
	end
	plot(g8[:,1], g8[:,3:end],ls=:dash, shape=:circle,markersize=2)
	plot!(title= "N=8", xlabel = "J/Γ", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ e2fee3ab-624f-4c97-83ad-921bf754483e
begin
	g10 = readdlm("./data/exact/exact-fixg10.txt")
	for i= 3:129
		g10[:,i] = g10[:,i] .- g10[:,2]
	end
	plot(g10[:,1], g10[:,3:end],ls=:dash, shape=:circle,markersize=2)
	plot!(title= "N=10", xlabel = "J/Γ", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ e907a3fd-7ce7-4d93-be4b-f18210be0386
begin
	g16 = readdlm("./data/exact/exact-fixg16.txt")
	for i= 3:129
		g16[:,i] = g16[:,i] .- g16[:,2]
	end
	plot(g16[:,1], g16[:,3:end],ls=:dash, shape=:circle,markersize=2)
	plot!(title= "N=16", xlabel = "J/Γ", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ 0ee9d481-6b5f-4ff0-971e-b3377af42179
begin
	j7 = readdlm("./data/exact/exact-fixj7.txt")
	for i= 3:129
		j7[:,i] = j7[:,i] .- j7[:,2]
	end
	plot(j7[:,1], j7[:,3:end],ls=:dash, shape=:circle,markersize=2)
	plot!(title= "N=7", xlabel = "Γ/J", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ 7f1a793e-ac02-4e35-b692-5df7681fef83
begin
	j8 = readdlm("./data/exact/exact-fixj8.txt")
	for i= 3:129
		j8[:,i] = j8[:,i] .- j8[:,2]
	end
	plot(j8[:,1], j8[:,3:end],ls=:dash, shape=:circle,markersize=2)
	plot!(title= "N=8", xlabel = "Γ/J", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ cbd44adc-2b17-49c0-b9a3-4ad0a6b57504
begin
	j10 = readdlm("./data/exact/exact-fixj10.txt")
	for i= 3:129
		j10[:,i] = j10[:,i] .- j10[:,2]
	end
	plot(j10[:,1], j10[:,3:end],ls=:dash, shape=:circle,markersize=2)
	plot!(title= "N=10", xlabel = "Γ/J", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ 4a2778a3-eeb2-43a3-af07-09e75b3334e7
begin
	j16 = readdlm("./data/exact/exact-fixj16.txt")
	for i= 3:129
		j16[:,i] = j16[:,i] .- j16[:,2]
	end
	plot(j16[:,1], j16[:,3:end],ls=:dash, shape=:circle,markersize=2)
	plot!(title= "N=16", xlabel = "Γ/J", ylabel = "ΔEm",legend=:none)
end

# ╔═╡ Cell order:
# ╟─27c49980-570a-4845-9d7a-7c540f0d2abf
# ╟─166f0878-1ed4-4ca1-8118-fb650d65acce
# ╟─e2fee3ab-624f-4c97-83ad-921bf754483e
# ╟─e907a3fd-7ce7-4d93-be4b-f18210be0386
# ╠═0ee9d481-6b5f-4ff0-971e-b3377af42179
# ╠═7f1a793e-ac02-4e35-b692-5df7681fef83
# ╠═cbd44adc-2b17-49c0-b9a3-4ad0a6b57504
# ╠═4a2778a3-eeb2-43a3-af07-09e75b3334e7
# ╟─66817958-c4da-11eb-04ec-5f000da4c04c
