### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 9555258c-5890-11ec-1c67-1dba0d6df118
begin
	using Pkg; Pkg.activate("./")
	using cMPO
	using JLD, HDF5
	using Plots; gr(xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=12)
	default(palette = palette(:okabe_ito))
	using DelimitedFiles
	using Printf
	using PlutoUI
	"packages"
end

# ╔═╡ 3f602ef1-f8e8-4833-bd49-00c0e7caedc8
md"""
### Setups
"""

# ╔═╡ c9f8e371-5a63-4f84-a69f-6a2d0087d77c
invT = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0]

# ╔═╡ 7191ed9d-7be8-4f15-a262-51af41ff7180
lt = length(invT);

# ╔═╡ 79364cf8-a5b0-470f-9f93-7d700e2cb65f
N = 40

# ╔═╡ abf0db31-a739-4992-bd4f-10e81238aeef
D = 8

# ╔═╡ 670a4cd0-0db8-4f03-a24d-68aa5dc28bdf
begin
	pm = make_operator(pauli(:+),D)
	pz = make_operator(pauli(:z),D)
	sz = make_operator(0.5*pauli(:z),D)
	"pm, pz, sz"
end

# ╔═╡ be3eaa93-b4fc-4300-9455-00d4bf8be7e0
md"""
## Heisenberg model
"""

# ╔═╡ f5826335-5d41-4637-8dd8-5b35deda26ba
begin
	p_ReG = [@sprintf "./data/xxz/imagtime/giwn/Jz_1.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	p_ImG = [@sprintf "./data/xxz/imagtime/gdivwn/Jz_1.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	p_dReG = [@sprintf "./data/xxz/imagtime/dReG/Jz_1.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	
	d_ReG = [readdlm(p_ReG[i]) for i=1:lt]
	d_ImG = [readdlm(p_ImG[i]) for i=1:lt]
	d_dReG = [readdlm(p_dReG[i]) for i=1:lt]
	
	ReG = Vector{Matrix}(undef,lt)
	ImG_dE = Vector{Matrix}(undef,lt)
	dReG = Vector{Matrix}(undef,lt)
	
	for i = 1:lt
		b = invT[i]
	    ReG[i], ImG_dE[i], dReG[i] = zeros(N,2), zeros(N,2), zeros(N,2)
		ReG[i][:,1] = d_ReG[i][:,1]; ReG[i][:,2] = d_ReG[i][:,2]
		ImG_dE[i][:,1] = d_ImG[i][:,1]; ImG_dE[i][:,2] = -2/b*d_ImG[i][:,3]
		dReG[i][:,1] = d_dReG[i][:,1]; dReG[i][:,2] = 2/b*d_dReG[i][:,2]
	end
	"load ls data pm: ReG, ImG_dE, dReG"
end

# ╔═╡ 610a2f79-2762-4420-bf4f-fa53a6ad46dd
begin
	p_ReGz = [@sprintf "./data/xxz/imagtime/giwn/Jz_1.0_pz_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	p_ImGz = [@sprintf "./data/xxz/imagtime/gdivwn/Jz_1.0_pz_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	p_dReGz = [@sprintf "./data/xxz/imagtime/dReG/Jz_1.0_pz_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	
	d_ReGz = [readdlm(p_ReGz[i]) for i=1:lt]
	d_ImGz = [readdlm(p_ImGz[i]) for i=1:lt]
	d_dReGz = [readdlm(p_dReGz[i]) for i=1:lt]
	
	ReGz = Vector{Matrix}(undef,lt)
	ImG_dEz = Vector{Matrix}(undef,lt)
	dReGz = Vector{Matrix}(undef,lt)
	
	for i = 1:lt
		b = invT[i]
	    ReGz[i], ImG_dEz[i], dReGz[i] = zeros(N,2), zeros(N,2), zeros(N,2)
		ReGz[i][:,1] = d_ReGz[i][:,1]; ReGz[i][:,2] = d_ReGz[i][:,2]
		ImG_dEz[i][:,1] = d_ImGz[i][:,1]; ImG_dEz[i][:,2] = -2/b*d_ImGz[i][:,3]
		dReGz[i][:,1] = d_dReGz[i][:,1]; dReGz[i][:,2] = 2/b*d_dReGz[i][:,2]
	end
	"load ls data pz: ReG, ImG_dE, dReG"
end

# ╔═╡ 286b0bf7-278d-4799-bcd3-74a9bad0d5b4
begin
	data_p = @sprintf "./data/xxz/Jz_1.0_D_%i.jld" D 
	data = load(data_p)
	"load heisenberg cmpo data"
end

# ╔═╡ be07f683-dc13-4df5-9dde-c1d1e7d7acdb
@bind β Slider(invT)

# ╔═╡ aa1a850e-f3e3-4363-a8e6-fbe972289e15
t = findall(x->x == β, invT)[1];

# ╔═╡ 28550373-4cd4-4f8a-b1eb-1e3257d3ef0d
begin
	w = XXZmodel(1.0)
	key = string(β); ψ = tocmps(data[key][2])	
	Gt = [π*invT[i]/4*correlation_2time(invT[i]/2,pm,pm',ψ,w,invT[i]) for i=1:lt]
	Gtz = [π*invT[i]/4*correlation_2time(invT[i]/2,sz,sz',ψ,w,invT[i]) for i=1:lt]
	"πβ/4*G(β/2)"
end

# ╔═╡ 4c223588-5dff-45c2-bd2b-8910545946eb
begin
	scatter([0.0], [Gt[t]],label="πβ/4*G(β/2)")
	plot!(ImG_dE[t][:,1], ImG_dE[t][:,2], 
		line=(:dash),
		marker=(:circle),
		label="ImG(iωn)/ΔE")
	plot!(dReG[t][:,1],dReG[t][:,2], 
		line=(:dash),
		marker=(:circle),
		label="dReG_dωn")
	#plot!(xlim=(-0.5, xmax1))
	fig_heisenberg = @sprintf "AFM Heisenberg model, β=%i" β
	plot!(xlabel="ωn", ylabel="S(0)±",title=fig_heisenberg)
end

# ╔═╡ c41850cb-419c-4e27-87fd-a0bc9ad4644a
begin
	scatter([0.0], [Gtz[t]],label="πβ/4*G(β/2)")
	plot!(ImG_dEz[t][:,1], ImG_dEz[t][:,2]/4, 
		line=(:dash),
		marker=(:circle),
		label="ImG(iωn)/ΔE")
	plot!(dReGz[t][:,1],dReGz[t][:,2]/4, 
		line=(:dash),
		marker=(:circle),
		label="dReG_dωn")
	#plot!(xlim=(-0.5, xmax1))
	plot!(xlabel="ωn", ylabel="S(0)z",title=fig_heisenberg)
end

# ╔═╡ 25f6c7ed-6dc1-4a2c-8858-e1b9fbc940f3
@bind xmax1 Slider(ImG_dE[t][:,1])

# ╔═╡ 37cb3aff-f194-4b36-8bbc-bf5565d714cf
begin
	plot(ReG[t][:,1], ReG[t][:,2], 
		line=(:dash),
		marker=(:circle),
		label="ReG(iωn)")
	plot!(xlabel="ωn", ylabel="ReG(iωn)±",title=fig_heisenberg, legend=:bottomright)
end

# ╔═╡ 55b2e360-4f35-4ec6-94ad-65dcc789eb3f
begin
	plot(ReGz[t][:,1], ReGz[t][:,2], 
		line=(:dash),
		marker=(:circle),
		label="ReG(iωn)")
	plot!(xlabel="ωn", ylabel="ReG(iωn)z",title=fig_heisenberg, legend=:bottomright)
end

# ╔═╡ f1868b84-cb2a-478b-b348-d5ce0bdf3951
md"""
## XX model
"""

# ╔═╡ 59739cdc-027c-4526-a533-2211657087a2
begin
	p_ReG1 = [@sprintf "./data/xxz/imagtime/giwn/Jz_0.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	p_ImG1 = [@sprintf "./data/xxz/imagtime/gdivwn/Jz_0.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	p_dReG1 = [@sprintf "./data/xxz/imagtime/dReG/Jz_0.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	
	d_ReG1 = [readdlm(p_ReG1[i]) for i=1:lt]
	d_ImG1 = [readdlm(p_ImG1[i]) for i=1:lt]
	d_dReG1 = [readdlm(p_dReG1[i]) for i=1:lt]
	
	ReG1 = Vector{Matrix}(undef,lt)
	ImG_dE1 = Vector{Matrix}(undef,lt)
	dReG1 = Vector{Matrix}(undef,lt)
	
	for i = 1:lt
		b = invT[i]
	    ReG1[i], ImG_dE1[i], dReG1[i] = zeros(N,2), zeros(N,2), zeros(N,2)
		ReG1[i][:,1] = d_ReG1[i][:,1]; ReG1[i][:,2] = d_ReG1[i][:,2]
		ImG_dE1[i][:,1] = d_ImG1[i][:,1]; ImG_dE1[i][:,2] = -2/b*d_ImG1[i][:,3]
		dReG1[i][:,1] = d_dReG1[i][:,1]; dReG1[i][:,2] = 2/b*d_dReG1[i][:,2]
	end
	"load ls data pm: ReG, ImG_dE, dReG"
end

# ╔═╡ 8b8fc553-3ea7-4a69-9c2c-08c9b74753b3
begin
	p_ReG1z = [@sprintf "./data/xxz/imagtime/giwn/Jz_0.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	p_ImG1z = [@sprintf "./data/xxz/imagtime/gdivwn/Jz_0.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	p_dReG1z = [@sprintf "./data/xxz/imagtime/dReG/Jz_0.0_pm_D_%i_beta_%i.txt" D invT[i] for i=1:lt]
	
	d_ReG1z = [readdlm(p_ReG1z[i]) for i=1:lt]
	d_ImG1z = [readdlm(p_ImG1z[i]) for i=1:lt]
	d_dReG1z = [readdlm(p_dReG1z[i]) for i=1:lt]
	
	ReG1z = Vector{Matrix}(undef,lt)
	ImG_dE1z = Vector{Matrix}(undef,lt)
	dReG1z = Vector{Matrix}(undef,lt)
	
	for i = 1:lt
		b = invT[i]
	    ReG1z[i], ImG_dE1z[i], dReG1z[i] = zeros(N,2), zeros(N,2), zeros(N,2)
		ReG1z[i][:,1] = d_ReG1z[i][:,1]; ReG1z[i][:,2] = d_ReG1z[i][:,2]
		ImG_dE1z[i][:,1] = d_ImG1z[i][:,1]; ImG_dE1z[i][:,2] = -2/b*d_ImG1z[i][:,3]
		dReG1z[i][:,1] = d_dReG1z[i][:,1]; dReG1z[i][:,2] = 2/b*d_dReG1z[i][:,2]
	end
	"load ls data pz: ReG, ImG_dE, dReG"
end

# ╔═╡ e30a2f34-fc0c-4a0a-b6c9-9132c8d2ded3
begin
	data_p1 = @sprintf "./data/xxz/Jz_0.0_D_%i.jld" D 
	data1 = load(data_p1)
	"load heisenberg cmpo data"
end

# ╔═╡ 427ad555-3e33-4cd2-9a40-64e61254b73e
@bind β1 Slider(invT)

# ╔═╡ 55cdf228-e821-4620-a258-4e85f6ff693c
t1 = findall(x->x == β1, invT)[1];

# ╔═╡ c2991a77-8be4-4c49-84d8-63b7fe4ed897
begin
	w1 = XXZmodel(0.0)
	key1 = string(β1); ψ1 = tocmps(data1[key1][2])	
	Gt1 = [π*invT[i]/4*correlation_2time(invT[i]/2,pm,pm',ψ1,w1,invT[i]) for i=1:lt]
	Gt1z = [π*invT[i]/4*correlation_2time(invT[i]/2,sz,sz',ψ1,w1,invT[i]) for i=1:lt]
	"πβ/4*G(β/2)"
end

# ╔═╡ 0a55e47a-3048-494d-af2c-46a9ab9a0746
begin
	scatter([0.0], [Gt1[t1]],label="πβ/4*G(β/2)")
	plot!(ImG_dE1[t1][:,1], ImG_dE1[t1][:,2], 
		line=(:dash),
		marker=(:circle),
		label="ImG(iωn)/ΔE")
	plot!(dReG1[t1][:,1],dReG1[t1][:,2], 
		line=(:dash),
		marker=(:circle),
		label="dReG_dωn")
	#plot!(xlim=(-0.5, xmax1))
	fig_xx = @sprintf "XX model, β=%i" β1
	plot!(xlabel="ωn", ylabel="S(0)±",title=fig_xx)
end

# ╔═╡ 0d183ad6-5d76-4574-96e0-f35dfcafe6e6
begin
	scatter([0.0], [Gt1z[t1]],label="πβ/4*G(β/2)")
	plot!(ImG_dE1z[t1][:,1], ImG_dE1z[t1][:,2]/4, 
		line=(:dash),
		marker=(:circle),
		label="ImG(iωn)/ΔE")
	plot!(dReG1z[t1][:,1],dReG1z[t1][:,2]/4, 
		line=(:dash),
		marker=(:circle),
		label="dReG_dωn")
	#plot!(xlim=(-0.5, xmax1))
	plot!(xlabel="ωn", ylabel="S(0)z",title=fig_xx)
end

# ╔═╡ 45a20b84-e958-4d25-aab6-006abd54db18
begin
	plot(ReG1[t1][:,1], ReG1[t1][:,2], 
		line=(:dash),
		marker=(:circle),
		label="ReG(iωn)")
	plot!(xlabel="ωn", ylabel="ReG(iωn)±",title=fig_xx, legend=:bottomright)
end

# ╔═╡ 29785e0f-92ca-4085-bb34-1e4c79f5391e
begin
	plot(ReG1z[t1][:,1], ReG1z[t1][:,2], 
		line=(:dash),
		marker=(:circle),
		label="ReG(iωn)")
	plot!(xlabel="ωn", ylabel="ReG(iωn)z",title=fig_xx, legend=:bottomright)
end

# ╔═╡ Cell order:
# ╟─7191ed9d-7be8-4f15-a262-51af41ff7180
# ╟─aa1a850e-f3e3-4363-a8e6-fbe972289e15
# ╠═55cdf228-e821-4620-a258-4e85f6ff693c
# ╟─3f602ef1-f8e8-4833-bd49-00c0e7caedc8
# ╠═c9f8e371-5a63-4f84-a69f-6a2d0087d77c
# ╟─79364cf8-a5b0-470f-9f93-7d700e2cb65f
# ╟─abf0db31-a739-4992-bd4f-10e81238aeef
# ╟─670a4cd0-0db8-4f03-a24d-68aa5dc28bdf
# ╟─be3eaa93-b4fc-4300-9455-00d4bf8be7e0
# ╟─f5826335-5d41-4637-8dd8-5b35deda26ba
# ╟─610a2f79-2762-4420-bf4f-fa53a6ad46dd
# ╟─286b0bf7-278d-4799-bcd3-74a9bad0d5b4
# ╟─28550373-4cd4-4f8a-b1eb-1e3257d3ef0d
# ╠═be07f683-dc13-4df5-9dde-c1d1e7d7acdb
# ╟─c41850cb-419c-4e27-87fd-a0bc9ad4644a
# ╟─4c223588-5dff-45c2-bd2b-8910545946eb
# ╟─25f6c7ed-6dc1-4a2c-8858-e1b9fbc940f3
# ╟─37cb3aff-f194-4b36-8bbc-bf5565d714cf
# ╟─55b2e360-4f35-4ec6-94ad-65dcc789eb3f
# ╟─f1868b84-cb2a-478b-b348-d5ce0bdf3951
# ╟─59739cdc-027c-4526-a533-2211657087a2
# ╟─8b8fc553-3ea7-4a69-9c2c-08c9b74753b3
# ╟─e30a2f34-fc0c-4a0a-b6c9-9132c8d2ded3
# ╟─c2991a77-8be4-4c49-84d8-63b7fe4ed897
# ╠═427ad555-3e33-4cd2-9a40-64e61254b73e
# ╟─0a55e47a-3048-494d-af2c-46a9ab9a0746
# ╟─0d183ad6-5d76-4574-96e0-f35dfcafe6e6
# ╟─45a20b84-e958-4d25-aab6-006abd54db18
# ╟─29785e0f-92ca-4085-bb34-1e4c79f5391e
# ╟─9555258c-5890-11ec-1c67-1dba0d6df118
