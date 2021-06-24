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
function plot_exact(N::Int, xlab::Char; Nstates::Int = 128)
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
	Nstates += 1
	for i= 3:Nstates
		d[:,i] = d[:,i] .- d[:,2]
	end
	#plot(d[:,1], d[:,3:end] .* N ,ls=:dash, shape=:circle,markersize=2)
	plot(d[:,1], d[:,3] .* N, line = (:black, 1), 
		label ="ED")
	plot!(d[:,1], d[:,4:Nstates] .* N ,line = (:black, 1), 
		label = false, fmt = :png)
	t = @sprintf "N=%i" N
	plot!(title= t, xlabel = xl, ylabel = "ΔEm")	
end

# ╔═╡ e18f2569-480a-44f7-be42-046f28ce212c
function add_plot_analytic(N::Int, xlab::Char; Nstates::Int = 128)
	if xlab == 'Γ'
		path = @sprintf "../data/exact/analytic-fixj%i.txt" N
		xl = @sprintf "Γ/J"
	elseif xlab == 'J'
		path = @sprintf "../data/exact/analytic-fixg%i.txt" N
		xl = @sprintf "J/Γ"
	else 
		@error " xlabel should be 'Γ' or 'J'. "
	end
	Nstates += 1
	d = readdlm(path)
	for i= 3:Nstates
		d[:,i] = d[:,i] .- d[:,2]
	end
	#plot(d[:,1], d[:,3:end] .* N ,ls=:dash, shape=:circle,markersize=2)
	plot!(d[:,1], d[:,3], line = (:dash),
		 marker = (:circle, 2), label ="Analytic")
	plot!(d[:,1], d[:,4:Nstates], line = (:dash),
		 marker = (:circle, 2), label =false)
end

# ╔═╡ 27c49980-570a-4845-9d7a-7c540f0d2abf
begin
	plot_exact(8, 'J', Nstates = 128 )
	add_plot_analytic(8, 'J', Nstates = 110)
	plot!(legend = :topleft)
end

# ╔═╡ e907a3fd-7ce7-4d93-be4b-f18210be0386
begin
	plot_exact(16, 'J', Nstates = 128)
	add_plot_analytic(16, 'J', Nstates = 128)
	plot!(legend = :topright)
end

# ╔═╡ a1175c8f-ce31-456f-aefb-77360f5356f1
md"""
### Ising spectrum
"""

# ╔═╡ Cell order:
# ╠═27c49980-570a-4845-9d7a-7c540f0d2abf
# ╟─e907a3fd-7ce7-4d93-be4b-f18210be0386
# ╟─aae8ae23-b760-45cb-8541-8f3e650a97b3
# ╠═e18f2569-480a-44f7-be42-046f28ce212c
# ╟─66817958-c4da-11eb-04ec-5f000da4c04c
# ╟─a1175c8f-ce31-456f-aefb-77360f5356f1
# ╠═9e28f94d-72ae-4838-8ced-b2fd93ff3292
