### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ b77c1e29-028d-4a96-aa0f-a33644380b84
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ ce3a64cd-f495-486d-a67f-b5f318e746ad
begin
	using cMPO
	using DelimitedFiles, HDF5, JLD
	using Plots
	using Printf
	using LsqFit
	using LinearAlgebra
	using Arpack
end

# ╔═╡ 2aca9ac6-48ef-470f-a2ce-873b1acd7438
md"""
## cMPO results of transverse field Ising model
Setups: $\chi = 8$

gamma range: $[0.1, 0.5, 0.9, 1, 1.1, 2, 4]$

beta range: $[0.1, 1, 10, 20]$
"""

# ╔═╡ daeac86c-674b-411e-a354-94d30c2ccfb3
md"""
### Compare to exact results: f and $\langle \sigma_i^x \rangle$
"""

# ╔═╡ b2678f29-f16e-400f-a8d3-1c67dd5c258b
begin
	β = 20
	p1 = @sprintf "./data/f_and_sx_bnew_%.1f_s.txt" β
	d1 = readdlm(p1)
	p1
end

# ╔═╡ 65e5f4df-887e-4d32-be0b-0ffac4fd0651
begin
	t1_f = @sprintf "free energy loss at β=%.1f" β
	plot(d1[:,1], d1[:,3]-d1[:,2], ls=:dash, shape=:circle, label="floss")
	plot!(title= t1_f, xlabel = "Γ", ylabel = "floss",legend=:best)
end

# ╔═╡ a69d2c85-65dc-486f-8006-bde6720a648a
begin
	t1_x = @sprintf "Thermal average <σx> at β=%.1f" β
	plot(d1[:,1], d1[:,5], ls=:dash, shape=:circle, label="cMPO")
	plot!(d1[:,1], d1[:,4], lw=3,label="exact")
	plot!(title= t1_x, xlabel = "Γ", ylabel = "< σx >",legend=:bottomright)
end

# ╔═╡ 6155dce8-2737-44c9-b5d8-e2adeb905e6b
md"""
### Phase transition: $\langle \sigma_i^z \rangle$
Applied $z$-field $= 0.01$
"""

# ╔═╡ 3b0fe4ce-53f4-4014-9ed6-ea59dcad5d96
begin
	dz = readdlm("./data/sz_bnew_20.0_s2.txt")
	"./data/sz_bnew_20.0_s2.txt"
end

# ╔═╡ fdb7da5b-e32e-4f2a-9717-982faf99eaec
begin
	t1_z = @sprintf "Thermal average <σz> at β=%.1f" β
	plot(dz[:,1], dz[:,2], ls=:dash, shape=:circle, label="cMPO")
	plot!(title= t1_z, xlabel = "Γ", ylabel = "< σz >")
end

# ╔═╡ d5e7d64d-251d-4c56-b05b-d916e33b8c15
md"""
### Quantum critical point $\Gamma/J = 1$
"""

# ╔═╡ b66b22d1-f1fd-418a-95c3-69ce7e4b4847
Wc = TFIsing(1.0,1.0)

# ╔═╡ cac14c6e-3f0c-4726-91dc-244c52f10bee
#load data
begin
	sd = load("./data/gnew_1.0.jld")
	"./data/gnew_1.0.jld"
end

# ╔═╡ 6a715819-ae33-4868-bd9e-ddcc5f2d3894
beta = [i for i in range(0.1,30,step = 0.5)]

# ╔═╡ 0b6ce7d2-9efb-4248-ae7b-b5842d2eeab9
md"""
### NMR 1/T1
"""

# ╔═╡ 25df64ee-40aa-423d-b533-e19e912791a6
#load data
begin
	s1 = load("./data/tnew_0.5.jld")
	"./data/tnew_0.5.jld"
end

# ╔═╡ ceefb4a3-b1c3-4ff6-8970-a06e2f3e214f
#load data
begin
	s2 = load("./data/tnew_2.0.jld")
	"./data/tnew_2.0.jld"
end

# ╔═╡ 5380b277-cad3-49ca-b685-7e7060176ea4
md"""
### Informations in the specturm of $\vdash + \dashv$
Gap funtion $\Delta(\Gamma/J)$ at low temperature
"""

# ╔═╡ 9d2de768-a002-4ac9-949e-aa3b34d96be5
begin
	p2 = @sprintf "./data/bnew_%.1f_s.jld" β
	d2 = load(p2)
	p2
end

# ╔═╡ 2049cdf6-33aa-4b45-a910-a96e29629ad3
gamma = [i for i in range(0,3,step = 0.02)]

# ╔═╡ d4e3e159-2b83-4f02-bc9f-472b3dbc380e
begin
	p3 = @sprintf "./data/jnew_%.1f.jld" β
	d3 = load(p3)
	p3
end

# ╔═╡ 3e68755f-097f-4f4b-8cda-0c94f206ea98
J = [i for i in range(0.,1.,step = 0.02)]

# ╔═╡ d0498e1b-8987-4651-bba6-84f4002e2ab3
md"""
#### Partitian $\tilde{Z}'$
"""

# ╔═╡ c79da5e4-fdee-460a-ba69-36c546a4d15a
begin # Partitian function Z'
	Zc = zeros(length(beta))
	for i = 1:length(beta)
		b = beta[i]; k = string(beta[i])
		ψt = tocmps(sd[k][2])
		Zc[i] = partitian!(ψt,Wc,b).res
	end
	#plot
	t1_Zc = @sprintf "Z' at Γ/J = 1"
	plot(1 ./beta[5:end], Zc[5:end], ls=:dash, shape=:circle, label="Γ/J = 1.0")
	plot!(title= t1_Zc, xlabel = "T", ylabel = "Z'",legend=:topleft)
end

# ╔═╡ 08e629d6-05d3-4ae6-8def-0c479957244e
# number of eigen-states
num = 10

# ╔═╡ 8949d60a-0a6a-4821-8b53-0217b146ddf0
begin
	z = pauli('z')
	eye = Matrix(1.0I, 8, 8)
	sz = eye ⊗ z ⊗ eye
	"define pauli z"
end

# ╔═╡ 23f215ab-ab5f-4523-b5d5-44b8630c4eaf
begin
	βc = beta[41]; kc = string(βc) #select β
	ψc = tocmps(sd[kc][2])
	#cmpo
	tau_num = [i for i in range(0,βc,length=50)]
	Cz_num = [correlation_2time(x,z,z,ψc, Wc, βc) for x in tau_num]
	#theoretical
	tau_th = [i for i in range(0.05,βc-0.05,length=50)]
	Cz_th = [critical_zz_cor(x,βc) for x in tau_th]
	#plot
	t1_cz = @sprintf "< σz(τ,0) σz(0,0) > at β=%.1f" βc
	plot(tau_th, Cz_th, lw=2,label="theoretical")
	plot!(tau_num, Cz_num, ls=:dash, shape=:circle, label="cMPO")
	plot!(title= t1_cz, xlabel = "τ", ylabel = "< σz(τ,0) σz(0,0) >",legend=:best)
end

# ╔═╡ 90e412a0-9548-4270-8771-367546334c7f
begin
		#cmpo
		ω = [10^i for i in range(-3,1,length=100)]
		X_num1 = [imag_susceptibility(x,z,z,ψc, Wc, βc, η=0.05) for x in ω]
		X_num2 = [imag_susceptibility(x,z,z,ψc, Wc, βc, η=0.1) for x in ω]
		#theoretical
		X_th = [critical_zz_chi(x,βc) for x in ω]
		#plot
		t1_X = @sprintf "Im Χ(ω) at β=%.1f" βc
		plot(ω, X_th, lw=2,label="theoretical")
		plot!(ω, X_num1, ls=:dash, shape=:circle, label="cMPO:η=0.05")
		plot!(ω, X_num2, ls=:dash, shape=:circle, label="cMPO:η=0.1")
		plot!(title= t1_X, xlabel = "ω", ylabel = "Im Χ(ω)",legend=:best,
		xscale=:log)
end

# ╔═╡ 9004033e-e0be-4cfe-8f82-bb2d744865e7
begin #NMR 1/T1
	T1 = zeros(length(beta),2)
	for i = 1:length(beta)
		b = beta[i]; k = string(beta[i])
		ψt = tocmps(sd[k][2])
		T1[i,1] = structure_factor(0,z,z,ψt, Wc, b, η=0.05)
		T1[i,2] = structure_factor(0,z,z,ψt, Wc, b, η=0.1) 
	end
	T1_th = [structure_factor(1,1,x) for x in beta]
	#plot
	t1_T = @sprintf "NMR 1/T1 = S(0,0) at Γ/J = 1"
	plot(1 ./ beta, T1[:,1], ls=:dash, shape=:circle, label="cMPO:η=0.05")
	plot!(1 ./ beta, T1[:,2], ls=:dash, shape=:circle, label="cMPO:η=0.1")
	plot!(1 ./ beta, T1_th, lw=2,label="theoretical")
	plot!(title= t1_T, xlabel = "T", ylabel = "NMR 1/T1",legend=:best, xlim=(0,1))
end

# ╔═╡ 7db823ce-b609-43d6-a83b-42b99eec0566
begin #NMR 1/T1
	temp=[i for i in range(3,0.05,length=60)]
	T2 = zeros(length(temp),2)
	W1 = TFIsing(1.0,0.5)
	for i = 1:length(temp)
		b = 1/temp[i]; k = string(temp[i])
		ψt = tocmps(s1[k][2])
		T2[i,1] = structure_factor(0,z,z,ψt, W1, b, η=0.05)
		T2[i,2] = structure_factor(0,z,z,ψt, W1, b, η=0.1) 
	end
	#plot
	t1_T2 = @sprintf "NMR 1/T1 = S(0,0) at Γ/J = 0.5"
	plot(temp, T2[:,1], ls=:dash, shape=:circle, label="cMPO:η=0.05")
	plot!(temp, T2[:,2], ls=:dash, shape=:circle, label="cMPO:η=0.1")
	plot!(title= t1_T2, xlabel = "T", ylabel = "NMR 1/T1",legend=:best)
end

# ╔═╡ 43deae38-10b5-43ea-ad7b-8931b9edd02a
begin #partitian Z' 
	Z1 = zeros(length(temp))
	for i = 1:length(temp)
		b = 1/temp[i]; k = string(temp[i])
		ψt = tocmps(s1[k][2])
		Z1[i] = partitian!(ψt,W1,b).res
	end
	#plot
	t1_Z1 = @sprintf "Partitian Z' at Γ/J = 0.5"
	plot(temp[41:end], Z1[41:end], ls=:dash, shape=:circle, label="Γ/J = 2.0")
	#plot(temp, T3[:,2], ls=:dash, shape=:circle, label="cMPO:η=0.1")
	plot!(title= t1_Z1, xlabel = "T", ylabel = "Z'",legend=:topleft)
end

# ╔═╡ f0857adc-5155-4c9e-93df-16963b56413d
begin #NMR 1/T1
	T3 = zeros(length(temp),2)
	W2 = TFIsing(1.0,2.0)
	for i = 1:length(temp)
		b = 1/temp[i]; k = string(temp[i])
		ψt = tocmps(s2[k][2])
		T3[i,1] = structure_factor(0,z,z,ψt, W1, b, η=0.05)
		T3[i,2] = structure_factor(0,z,z,ψt, W1, b, η=0.1) 
	end
	#plot
	t1_T3 = @sprintf "NMR 1/T1 = S(0,0) at Γ/J = 2.0"
	plot(temp, T3[:,1], ls=:dash, shape=:circle, label="cMPO:η=0.05")
	#plot(temp, T3[:,2], ls=:dash, shape=:circle, label="cMPO:η=0.1")
	plot!(title= t1_T3, xlabel = "T", ylabel = "NMR 1/T1",legend=:topleft)
end

# ╔═╡ 5a453a6a-1c5e-4456-88f9-2bceec402be3
begin #partitian Z' 
	Z2 = zeros(length(temp))
	for i = 1:length(temp)
		b = 1/temp[i]; k = string(temp[i])
		ψt = tocmps(s2[k][2])
		Z2[i] = partitian!(ψt,W2,b).res
	end
	#plot
	t1_Z2 = @sprintf "Partitian Z' at Γ/J = 2.0"
	plot(temp[41:end], Z2[41:end], ls=:dash, shape=:circle, label="Γ/J = 2.0")
	#plot(temp, T3[:,2], ls=:dash, shape=:circle, label="cMPO:η=0.1")
	plot!(title= t1_Z2, xlabel = "T", ylabel = "Z'",legend=:topleft)
end

# ╔═╡ b98e458f-5b48-41f3-9b9c-f19335024513
begin
	lg = length(gamma)
	eng = zeros(lg,num) # eigen-energy
	ss = zeros(lg,num,num) # |Sz|
	for i = 1:lg
		g = gamma[i]; kg = string(g)
		ψg = tocmps(d2[kg][2])
		Wg = TFIsing(1.0,g)
		Kg = ψg * Wg * ψg
		e, v = eigen(Kg)
		s = v' * sz * v
		ss[i,:,:] = s[1:num,1:num]
		eng[i,:] = e[1:num] .-e[1]
	end	
end

# ╔═╡ d247d358-4367-4326-bcda-a11f0e66b83d
begin #plot energy level
	t1_e = @sprintf "Gap function at β=%.1f" β
	plot(gamma, eng, ls=:dash, shape=:circle,markersize=2)
	#plot(temp, T3[:,2], ls=:dash, shape=:circle, label="cMPO:η=0.1")
	plot!(gamma, 2 .* abs.(gamma .- 1.0),lw = 2, label=("|Δ|"))
	plot!(gamma, 4 .* abs.(gamma .- 1.0),lw = 2, label=("2|Δ|"))
	plot!(title= t1_e, xlabel = "Γ/J", ylabel = "ΔEm",legend=:topleft)
end

# ╔═╡ 1659f997-95b4-4a25-b854-ff33def73931
begin
	lj = length(J)
	eng2 = zeros(lj,num) # eigen-energy
	ss2 = zeros(lj,num,num) # |Sz|
	for i = 1:lj
		j = J[i]; kg = string(j)
		ψg = tocmps(d3[kg][2])
		Wg = TFIsing(j,1.0)
		Kg = ψg * Wg * ψg
		e, v = eigen(Kg)
		s = v' * sz * v
		ss2[i,:,:] = s[1:num,1:num]
		eng2[i,:] = e[1:num] .-e[1]
	end	
end

# ╔═╡ b69aa139-3ad5-4dc5-981a-bf0818656f21
begin #plot energy level
	t1_ej = @sprintf "Gap function in paramagnetic phase at β=%.1f" β
	plot(J, eng2, ls=:dash, shape=:circle,markersize=2)
	#plot(temp, T3[:,2], ls=:dash, shape=:circle, label="cMPO:η=0.1")
	plot!(J, 2 .* abs.(J .- 1.0),lw = 2, label=("|Δ|"))
	plot!(J, 4 .* abs.(J .- 1.0),lw = 2, label=("2|Δ|"))
	plot!(title= t1_ej, xlabel = "J/Γ", ylabel = "ΔEm",legend=:topright)
end

# ╔═╡ cace25bf-b722-4982-92c6-12095d1aecc2
begin
	lb = length(beta)
	engc = zeros(lb,num) # eigen-energy
	ssc = zeros(lb,num,num) # |Sz|
	for i = 1:lb
		b = beta[i]; kg = string(b)
		ψg = tocmps(sd[kg][2])
		Kg = ψg * Wc * ψg
		e, v = eigen(Kg, sortby=x->real(x))
		s = v' * sz * v
		ssc[i,:,:] = s[1:num,1:num]
		engc[i,:] = (e[1:num] .-e[1]) * b
	end	
end

# ╔═╡ ad6d3ec4-123d-453f-aead-c7072b11de4d
begin #plot energy level
	t1_ec = @sprintf "ΔE/T ~ -ln(T)"
	px = log.(beta)
	plot(px, engc, ls=:dash, shape=:circle,markersize=2)
	plot!(px, 0.25 .* px .+ 0.7,lw = 1.5, c=:blue,label=("-0.25log(T)"))
	plot!(title= t1_ec, xlabel = "-ln(T)", ylabel = "ΔEm/T",legend=:topleft)
end

# ╔═╡ d25175e8-b918-11eb-2ca0-2168adee0ebe
pwd()

# ╔═╡ Cell order:
# ╟─2aca9ac6-48ef-470f-a2ce-873b1acd7438
# ╟─daeac86c-674b-411e-a354-94d30c2ccfb3
# ╟─b2678f29-f16e-400f-a8d3-1c67dd5c258b
# ╠═65e5f4df-887e-4d32-be0b-0ffac4fd0651
# ╠═a69d2c85-65dc-486f-8006-bde6720a648a
# ╟─6155dce8-2737-44c9-b5d8-e2adeb905e6b
# ╟─3b0fe4ce-53f4-4014-9ed6-ea59dcad5d96
# ╟─fdb7da5b-e32e-4f2a-9717-982faf99eaec
# ╟─d5e7d64d-251d-4c56-b05b-d916e33b8c15
# ╠═b66b22d1-f1fd-418a-95c3-69ce7e4b4847
# ╠═cac14c6e-3f0c-4726-91dc-244c52f10bee
# ╠═6a715819-ae33-4868-bd9e-ddcc5f2d3894
# ╠═23f215ab-ab5f-4523-b5d5-44b8630c4eaf
# ╟─90e412a0-9548-4270-8771-367546334c7f
# ╟─9004033e-e0be-4cfe-8f82-bb2d744865e7
# ╟─0b6ce7d2-9efb-4248-ae7b-b5842d2eeab9
# ╟─25df64ee-40aa-423d-b533-e19e912791a6
# ╟─7db823ce-b609-43d6-a83b-42b99eec0566
# ╟─ceefb4a3-b1c3-4ff6-8970-a06e2f3e214f
# ╠═f0857adc-5155-4c9e-93df-16963b56413d
# ╟─5380b277-cad3-49ca-b685-7e7060176ea4
# ╟─9d2de768-a002-4ac9-949e-aa3b34d96be5
# ╠═2049cdf6-33aa-4b45-a910-a96e29629ad3
# ╠═b98e458f-5b48-41f3-9b9c-f19335024513
# ╠═d247d358-4367-4326-bcda-a11f0e66b83d
# ╠═d4e3e159-2b83-4f02-bc9f-472b3dbc380e
# ╠═3e68755f-097f-4f4b-8cda-0c94f206ea98
# ╠═1659f997-95b4-4a25-b854-ff33def73931
# ╠═b69aa139-3ad5-4dc5-981a-bf0818656f21
# ╟─5a453a6a-1c5e-4456-88f9-2bceec402be3
# ╟─43deae38-10b5-43ea-ad7b-8931b9edd02a
# ╟─d0498e1b-8987-4651-bba6-84f4002e2ab3
# ╠═c79da5e4-fdee-460a-ba69-36c546a4d15a
# ╠═cace25bf-b722-4982-92c6-12095d1aecc2
# ╠═ad6d3ec4-123d-453f-aead-c7072b11de4d
# ╠═08e629d6-05d3-4ae6-8def-0c479957244e
# ╟─8949d60a-0a6a-4821-8b53-0217b146ddf0
# ╠═ce3a64cd-f495-486d-a67f-b5f318e746ad
# ╟─b77c1e29-028d-4a96-aa0f-a33644380b84
# ╠═d25175e8-b918-11eb-2ca0-2168adee0ebe
