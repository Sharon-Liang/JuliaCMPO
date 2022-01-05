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
	using LinearAlgebra
	using LsqFit
	"packages"
end

# ╔═╡ be3eaa93-b4fc-4300-9455-00d4bf8be7e0
md"""
## Transverse filed Ising model
"""

# ╔═╡ 3f602ef1-f8e8-4833-bd49-00c0e7caedc8
md"""
### Setups
"""

# ╔═╡ c9f8e371-5a63-4f84-a69f-6a2d0087d77c
invT = [10.0, 20.0, 30.0, 40.0]

# ╔═╡ 7191ed9d-7be8-4f15-a262-51af41ff7180
lt = length(invT);

# ╔═╡ b719a415-ae92-4aa1-8597-6d373a76c608
gamma = [1.0]

# ╔═╡ 6ffb4dd0-0069-4924-8ec0-0782fcc96376
lg = length(gamma);

# ╔═╡ 79364cf8-a5b0-470f-9f93-7d700e2cb65f
N = 40

# ╔═╡ abf0db31-a739-4992-bd4f-10e81238aeef
D = 8

# ╔═╡ 42d2cb9d-96ba-46e6-88f0-74af89eb4511


# ╔═╡ 286b0bf7-278d-4799-bcd3-74a9bad0d5b4
begin
	data_path = [@sprintf "./data/ising/D_%i/g_%.1f.jld" D gamma[i] for i=1:lg]
	data = [load(data_path[i]) for i=1:lg]
	"load ising cmpo data"
end

# ╔═╡ 3d85febe-90aa-4713-af4d-1f3866941141
function diffaddexp(b::Real, e1::Real, e2::Real)
    if abs(e2 - e1) < 1.e-10
        return exp(-b*e1) * b
    else
        num = exp(-b*e1) - exp(-b*e2)
        den = e2 - e1
        return num/den
    end
end

# ╔═╡ 28d47fd0-11c7-47a8-9722-8e71d8c3e581
function Masubara_freq_GFdivOmega(z::Number, A::AbstractMatrix,B::AbstractMatrix,
                                ψ::CMPS, W::CMPO, β::Real)
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    #m = maximum(-β * e)
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    #den = exp.(-β * e .- m) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = A[i,j] * B[j,i] * diffaddexp(β,e[i],e[j])
        down = 1.0im * z - e[j] + e[i]
        num += up/down
    end
    return num/den
end

# ╔═╡ acd6afad-f3cc-4a55-8818-a7444a6a2b71
function ∂ReG(z::Number, A::AbstractMatrix,B::AbstractMatrix,
    ψ::CMPS, W::CMPO, β::Real)
    λ = 1.0
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(K)
    min = minimum(e); e = e .- min
    A = v' * A * v
    B = v' * B * v
    den = exp.(-β * e) |> sum
    num = 0.0
    for i = 1: length(e), j = 1: length(e)
        up = exp(-β*e[i]) - λ*exp(-β*e[j]) 
        up = up * A[i,j] * B[j,i] * (e[i] - e[j]) * (-2z)
        down = z^2 + (e[i] - e[j])^2
        down = down^2
        num += up/down
    end
    return num/den
end

# ╔═╡ 2e95d260-31fd-45ea-8da3-d7f734f4b590


# ╔═╡ dea97d84-956d-47d6-b563-bb2287b38585
@bind g Slider(gamma)

# ╔═╡ 297dbb2d-18b1-4fbb-be75-85a4ed256a21
gind = findall(x->x == g, gamma)[1];

# ╔═╡ 918911a0-a78c-4713-84f2-16d901b72adf
begin
	s0, s1 = zeros(length(invT)), zeros(length(invT))
	if g == 1.0
		pl0 = [@sprintf "./data/ising/spectrum-Li/Sw/g_%.1f_beta_%i.txt" g invT[i] for i=1:length(invT)]
		pl1 = [@sprintf "./data/ising/spectrum-Li/Sw/g_%.1f_beta_%i_eta_2pidbeta.txt" g invT[i] for i=1:length(invT)]
		dl0 = [readdlm(pl0[i]) for i=1:lt]
		dl1 = [readdlm(pl1[i]) for i=1:lt]
		s0 =[dl0[i][800,2] for i=1:lt]
		s1 =[dl1[i][800,2] for i=1:lt]
	end
	"load Zi-long Li data: S(0,η=0.001), S(0,η=2π/β) "
end

<<<<<<< HEAD
# ╔═╡ c8f38627-b5dc-4f18-ae54-beb2406734f0
s0

# ╔═╡ 945e2062-13aa-4ba2-9e33-f030f9cdeeb2
begin
	plot(dl0[t][:,1], dl0[t][:,2], label=@sprintf "η=0.001, β=%i" β)
	plot!(ylabel="S", xlabel="ω")
end

# ╔═╡ 4812d37b-8283-48a2-8f61-b81b63c99c62
begin
	plot(dl1[t][:,1], dl1[t][:,2], label=@sprintf "η=2π/β, β=%i" β)
	plot!(ylabel="S", xlabel="ω")
end

=======
>>>>>>> 903bd3d (update .gitignore)
# ╔═╡ 9bf1b230-951d-40d6-a759-db7740490511
begin
	if g == 1.0
		pw1 = [@sprintf "./data/ising/imagtime/S0iwn/g_%.1f_D_%i_beta_%i.txt" g D invT[i] for i=1:length(invT)]
		pw2 = [@sprintf "./data/ising/imagtime/S0iwn/g_%.1f_D_%im2_beta_%i.txt" g D invT[i] for i=1:length(invT)]
		dw1 = [readdlm(pw1[i]) for i=1:lt]
		dw2 = [readdlm(pw2[i]) for i=1:lt]
	end
	"load cmpo data: S(0,η=ωn) "
end

<<<<<<< HEAD
# ╔═╡ a31b1f2f-3102-4e7e-a5a0-6003a785a3d5
begin
	plot(xlabel="2π/β", ylabel="S(0,η)", title=(@sprintf "g=%.1f, β=%i" g β))
	scatter!(dw1[t][:,1], dw1[t][:,2], label="cmpo χ=8")
	#scatter!(dw2[t][:,1], dw2[t][:,2], label="cmpo χ=8×2")
	scatter!([2π/β],[s1[t]],marker=(:star, stroke(0.2)), label="η=2π/β")
	scatter!([0.0],[s0[t]],marker=(:star, stroke(0.2)), label="η=0.001")
end

=======
>>>>>>> 903bd3d (update .gitignore)
# ╔═╡ ef94f39d-0ab5-4259-8531-5279433e7996
begin
	plot(xlabel="2π/β", ylabel="S(0,η)", title=(@sprintf "g=%.1f" g))
	for i = lt:-1:1
		plot!(dw1[i][:,1], dw1[i][:,2], 
			line=(:dash,0.5),
			marker=(:circle, 3, stroke(0.2)), 
			label=@sprintf "β=%i" invT[i])
		scatter!([2π/invT[i]],[s1[i]],marker=(:star, stroke(0.2)), label=false)
		scatter!([0.0],[s0[i]],marker=(:star, stroke(0.2)), label="2π/β")
	end
	plot!()
	#scatter!(dw2[t][:,1], dw2[t][:,2], label="cmpo χ=8×2")
	#scatter!([2π/β],[s1[t]],marker=(:star, stroke(0.2)), label="η=0.001")
	#scatter!([0.0],[s0[t]],marker=(:star, stroke(0.2)), label="2π/β")
end

# ╔═╡ f5826335-5d41-4637-8dd8-5b35deda26ba
begin
	p_ReG = [@sprintf "./data/ising/imagtime/giwn/g_%.1f_D_%i_beta_%i.txt" g D invT[i] for i=1:lt]
	p_ImG = [@sprintf "./data/ising/imagtime/gdivwn/g_%.1f_D_%i_beta_%i.txt" g D invT[i] for i=1:lt]
	p_dReG = [@sprintf "./data/ising/imagtime/dReG/g_%.1f_D_%i_beta_%i.txt" g D invT[i] for i=1:lt]
	
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
	"load ls data: ReG, ImG_dE, dReG"
end

# ╔═╡ dbb09e85-141f-42dc-8fe9-7fd3a1ddde95
path = [@sprintf "./data/ising/imagtime/dReG/g_%.1f_D_%i_beta_%i.txt" g D invT[i] for i=1:length(invT)]

# ╔═╡ e0718b1f-be57-47b4-adbb-fb6f0f626c91
d1 = [readdlm(path[i]) for i=1:length(invT)];

# ╔═╡ 4e8df580-1561-4506-ad1a-582e431e5caf
begin
	plot(d1[1][:,1], d1[1][:,2], marker=(:circle,2),label= @sprintf "β=%i" invT[1])
	for i=2:length(invT)
		plot!(d1[i][:,1], d1[i][:,2]./invT[i], marker=(:circle,2),label= @sprintf "β=%i" invT[i])
	end
	#plot!(d1[length(invT)][:,1], d1[length(invT)][:,2], marker=(:circle,2), label= @sprintf "β=%i" invT[length(invT)])
	plot!(xlim=(-0.1,5))
end

# ╔═╡ be07f683-dc13-4df5-9dde-c1d1e7d7acdb
@bind β Slider(invT)

# ╔═╡ aa1a850e-f3e3-4363-a8e6-fbe972289e15
t = findall(x->x == β, invT)[1];

# ╔═╡ a31b1f2f-3102-4e7e-a5a0-6003a785a3d5
begin
	plot(xlabel="2π/β", ylabel="S(0,η)", title=(@sprintf "g=%.1f, β=%i" g β))
	scatter!(dw1[t][:,1], dw1[t][:,2], label="cmpo χ=8")
	#scatter!(dw2[t][:,1], dw2[t][:,2], label="cmpo χ=8×2")
	scatter!(dl1[t][:,1],dl1[t][:,2],marker=(:star, stroke(0.2)), label="η=2πn/β")
	scatter!([0.0],[s0[t]],marker=(:star, stroke(0.2)), label="η=0.001")
end

# ╔═╡ 78307144-6fa2-40f8-b18c-63fe3c63fa28
begin
	plot(xlabel="2π/β", ylabel="S(0,η) loss", title=(@sprintf "g=%.1f, β=%i" g β))
	plot!(dw1[t][:,1], abs.(dw1[t][:,2] .- dl1[t][:,2]),
		line=(:dash,1),
		marker=(:circle,4, stroke(0.2)),
		label=false)
end

# ╔═╡ 28550373-4cd4-4f8a-b1eb-1e3257d3ef0d
begin
	w = TFIsing(1.0, g)
	key = string(β); ψ = tocmps(data[gind][key][2])	
	pz = make_operator(pauli(:z),D)
	Gt = [π*invT[i]/4*correlation_2time(invT[i]/2,pz,pz',ψ,w,invT[i]) for i=1:lt]
	"πβ/4*G(β/2)"
end

# ╔═╡ 7041282e-d0ed-41e0-a4cd-5bdd70339a81
gw = [Masubara_freq_GFdivOmega(i, pz, pz', ψ, w, β) for i in range(1.e-3,2,length=100)]

# ╔═╡ 8fb177fe-2598-4a8f-969c-175d1a9623d9
gw1 = [∂ReG(i, pz, pz', ψ, w, β) for i in range(1.e-3,2,length=100)]

# ╔═╡ 6746901a-03f3-4af2-bde7-752060a350a8
begin
	omega = range(1.e-3,2,length=100)
	plot(omega, -imag.(gw), label="ReG/dE")
	plot!(omega, gw1, label="dReG")
end

# ╔═╡ 4c223588-5dff-45c2-bd2b-8910545946eb
begin
	scatter([2π/β], [dl1[t][1,2]],label="Li, η=2π/β")
	scatter!([0.0], [Gt[t]],label="πβ/4*G(β/2)")
	plot!(ImG_dE[t][:,1], ImG_dE[t][:,2], 
		line=(:dash),
		marker=(:circle,stroke(0.1)),
		label="ImG(iωn)/ΔE")
	plot!(dReG[t][:,1],dReG[t][:,2], 
		line=(:dash),
		marker=(:circle,stroke(0.1)),
		label="dReG_dωn")
	plot!(ReG[t][:,1],-ReG[t][:,2] ./ ReG[t][:,1]*2/β, 
		line=(:dash),
		marker=(:star, stroke(0.1)),
		label="ReG/ωn")
	#plot!(xlim=(-0.3, xmax1))
	fig_s0 = @sprintf "ising g=%.1f, β=%i" g β
	plot!(xlabel="ωn", ylabel="S(0)",title=fig_s0)
end

# ╔═╡ 25f6c7ed-6dc1-4a2c-8858-e1b9fbc940f3
@bind xmax1 Slider(ImG_dE[t][:,1])

# ╔═╡ e2d89ddf-0958-46ae-b179-be6cd0a5c963


# ╔═╡ 37cb3aff-f194-4b36-8bbc-bf5565d714cf
begin
	plot(ReG[t][:,1], ReG[t][:,2], 
		line=(:dash),
		marker=(:circle),
		label="ReG(iωn)")
	fig_ReG = @sprintf "ising g=%.1f, β=%i" g β
	plot!(xlabel="ωn", ylabel="ReG(iωn)",title=fig_ReG, legend=:bottomright)
end

# ╔═╡ ba566041-7dd2-4bcf-adc7-d1479ed8f629


# ╔═╡ 8a2b197c-dca1-41b5-acf0-cd6f74f669a3
@bind xmax2 Slider(dReG[t][:,1])

# ╔═╡ Cell order:
# ╠═7191ed9d-7be8-4f15-a262-51af41ff7180
# ╠═6ffb4dd0-0069-4924-8ec0-0782fcc96376
# ╠═aa1a850e-f3e3-4363-a8e6-fbe972289e15
# ╠═297dbb2d-18b1-4fbb-be75-85a4ed256a21
# ╟─be3eaa93-b4fc-4300-9455-00d4bf8be7e0
# ╟─3f602ef1-f8e8-4833-bd49-00c0e7caedc8
# ╟─c9f8e371-5a63-4f84-a69f-6a2d0087d77c
# ╟─b719a415-ae92-4aa1-8597-6d373a76c608
# ╟─79364cf8-a5b0-470f-9f93-7d700e2cb65f
# ╟─abf0db31-a739-4992-bd4f-10e81238aeef
# ╟─918911a0-a78c-4713-84f2-16d901b72adf
# ╟─9bf1b230-951d-40d6-a759-db7740490511
<<<<<<< HEAD
# ╠═c8f38627-b5dc-4f18-ae54-beb2406734f0
# ╠═be07f683-dc13-4df5-9dde-c1d1e7d7acdb
# ╠═a31b1f2f-3102-4e7e-a5a0-6003a785a3d5
# ╠═ef94f39d-0ab5-4259-8531-5279433e7996
# ╟─945e2062-13aa-4ba2-9e33-f030f9cdeeb2
# ╟─4812d37b-8283-48a2-8f61-b81b63c99c62
=======
# ╠═a31b1f2f-3102-4e7e-a5a0-6003a785a3d5
# ╟─78307144-6fa2-40f8-b18c-63fe3c63fa28
# ╟─ef94f39d-0ab5-4259-8531-5279433e7996
>>>>>>> 903bd3d (update .gitignore)
# ╟─f5826335-5d41-4637-8dd8-5b35deda26ba
# ╠═dbb09e85-141f-42dc-8fe9-7fd3a1ddde95
# ╠═e0718b1f-be57-47b4-adbb-fb6f0f626c91
# ╠═4e8df580-1561-4506-ad1a-582e431e5caf
# ╠═42d2cb9d-96ba-46e6-88f0-74af89eb4511
# ╠═286b0bf7-278d-4799-bcd3-74a9bad0d5b4
# ╟─28550373-4cd4-4f8a-b1eb-1e3257d3ef0d
# ╟─3d85febe-90aa-4713-af4d-1f3866941141
# ╟─28d47fd0-11c7-47a8-9722-8e71d8c3e581
# ╟─7041282e-d0ed-41e0-a4cd-5bdd70339a81
# ╟─acd6afad-f3cc-4a55-8818-a7444a6a2b71
# ╠═8fb177fe-2598-4a8f-969c-175d1a9623d9
# ╟─6746901a-03f3-4af2-bde7-752060a350a8
# ╠═2e95d260-31fd-45ea-8da3-d7f734f4b590
# ╠═dea97d84-956d-47d6-b563-bb2287b38585
# ╟─be07f683-dc13-4df5-9dde-c1d1e7d7acdb
# ╠═4c223588-5dff-45c2-bd2b-8910545946eb
# ╟─25f6c7ed-6dc1-4a2c-8858-e1b9fbc940f3
# ╠═e2d89ddf-0958-46ae-b179-be6cd0a5c963
# ╠═37cb3aff-f194-4b36-8bbc-bf5565d714cf
# ╠═ba566041-7dd2-4bcf-adc7-d1479ed8f629
# ╠═8a2b197c-dca1-41b5-acf0-cd6f74f669a3
# ╟─9555258c-5890-11ec-1c67-1dba0d6df118
