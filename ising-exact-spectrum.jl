### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 66817958-c4da-11eb-04ec-5f000da4c04c
begin
	using DelimitedFiles
	using JLD, HDF5
	using Plots
	using Printf
end

# ╔═╡ 9e28f94d-72ae-4838-8ced-b2fd93ff3292
using cMPO

# ╔═╡ aae8ae23-b760-45cb-8541-8f3e650a97b3
function plot_exact(N::Integer, xlab::Char)
	if xlab == 'Γ'
		path = @sprintf "./data/exact/exact-fixj%i.txt" N
		xl = @sprintf "Γ/J"
	elseif xlab == 'J'
		path = @sprintf "./data/exact/exact-fixg%i.txt" N
		xl = @sprintf "J/Γ"
	else 
		@error " xlabel should be 'Γ' or 'J'. "
	end
	d = readdlm(path)
	for i= 3:129
		d[:,i] = d[:,i] .- d[:,2]
	end
	#plot(d[:,1], d[:,3:end] .* N ,ls=:dash, shape=:circle,markersize=2)
	plot(d[:,1], d[:,3:end] .* N ,lw=1, lc=:black)
	t = @sprintf "N=%i" N
	plot!(title= t, xlabel = xl, ylabel = "ΔEm",legend=:none)	
end

# ╔═╡ 166f0878-1ed4-4ca1-8118-fb650d65acce
plot_exact(7, 'J')

# ╔═╡ 27c49980-570a-4845-9d7a-7c540f0d2abf
begin
	plot_exact(8, 'J')
	ad8 = readdlm("./data/exact/analytic-fixg8.txt")
	for i= 3:size(ad8)[2]
		ad8[:,i] = ad8[:,i] .- ad8[:,2]
	end
	plot!(ad8[:,1], ad8[:,3:end] ,ls=:dash, shape=:circle,markersize=2)
	#plot!(ylims = [8,10])
end

# ╔═╡ e2fee3ab-624f-4c97-83ad-921bf754483e
plot_exact(10, 'J')

# ╔═╡ e907a3fd-7ce7-4d93-be4b-f18210be0386
begin
	plot_exact(16, 'J')
	ad16 = readdlm("./data/exact/analytic-fixg16.txt")
	for i= 3:size(ad16)[2]
		ad16[:,i] = ad16[:,i] .- ad16[:,2]
	end
	plot!(ad16[:,1], ad16[:,3:end] ,ls=:dash, shape=:circle,markersize=2)
	#plot!(ylims = [0,8])
end

# ╔═╡ 0ee9d481-6b5f-4ff0-971e-b3377af42179
plot_exact(7,'Γ')

# ╔═╡ 7f1a793e-ac02-4e35-b692-5df7681fef83
plot_exact(8,'Γ')

# ╔═╡ cbd44adc-2b17-49c0-b9a3-4ad0a6b57504
plot_exact(10,'Γ')

# ╔═╡ 4a2778a3-eeb2-43a3-af07-09e75b3334e7
plot_exact(16, 'Γ')

# ╔═╡ a1175c8f-ce31-456f-aefb-77360f5356f1
md"""
### Ising spectrum
"""

# ╔═╡ dfe3a339-ff34-4328-9b10-83dc6b31cb29
begin
	L = 200; J = 1; Γ = 5
	k = [π/L * x for x in range(-L,L-1,step=1)]
	e = [energy_density(x,J,Γ) for x in k]
	plot(k,e)
end

# ╔═╡ b08fa981-0f05-4fa0-b0d2-718ca1b226d8
ad100 = readdlm("./data/exact/analytic-fixg100.txt")

# ╔═╡ Cell order:
# ╠═166f0878-1ed4-4ca1-8118-fb650d65acce
# ╠═27c49980-570a-4845-9d7a-7c540f0d2abf
# ╠═e2fee3ab-624f-4c97-83ad-921bf754483e
# ╠═e907a3fd-7ce7-4d93-be4b-f18210be0386
# ╟─0ee9d481-6b5f-4ff0-971e-b3377af42179
# ╟─7f1a793e-ac02-4e35-b692-5df7681fef83
# ╟─cbd44adc-2b17-49c0-b9a3-4ad0a6b57504
# ╟─4a2778a3-eeb2-43a3-af07-09e75b3334e7
# ╟─aae8ae23-b760-45cb-8541-8f3e650a97b3
# ╟─66817958-c4da-11eb-04ec-5f000da4c04c
# ╟─a1175c8f-ce31-456f-aefb-77360f5356f1
# ╠═9e28f94d-72ae-4838-8ced-b2fd93ff3292
# ╠═dfe3a339-ff34-4328-9b10-83dc6b31cb29
# ╠═b08fa981-0f05-4fa0-b0d2-718ca1b226d8
