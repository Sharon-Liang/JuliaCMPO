### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 86574854-9b9c-4fa8-940f-20446705fce6
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 91909ae8-cb1f-4653-9ba6-144bd4d01747
using LinearAlgebra

# ╔═╡ 76d940bc-09c4-422c-ab5a-50c8be7534e2
begin
	using cMPO
	using DelimitedFiles, HDF5, JLD
	using Plots
	using Printf
	using LsqFit
end

# ╔═╡ a63b4ad6-46e6-4635-8037-c24b20abf070
# 矩阵元随gamma变化
# 矩阵元分析

# ╔═╡ 181610df-4507-4cec-a7e0-2f919aec32bd
begin
	# gamma = [0.1, 0.5, 0.9, 1, 1.1, 2, 4]
	g = 1.0
	path2 = @sprintf "./data/gnew_%.1f.jld" g
	g1 = load(path2)
end

# ╔═╡ da8e7582-314f-401b-9406-a334c299c6b2
beta = [i for i in range(0.1,30,step = 0.5)]

# ╔═╡ ddaf571f-8bea-4c47-b8f7-1d7499de4719
T = 1 ./ beta

# ╔═╡ 0cd9aa6f-c170-4fd1-9284-e431b76cd821
w1 = TFIsing(1.0, g)

# ╔═╡ b6015a68-602b-4883-9cc3-9e57c8dffdd7
begin
	# beta = [0.1, 1, 10, 20]
	b = 20
	path1 = @sprintf "./data/bnew_%.1f.jld" b
	p1 = @sprintf "./data/f_and_sx_bnew_%.1f.txt" b
	p2 = @sprintf "./data/f_and_sx_bnew_%.1f_r.txt" b
	p3 = @sprintf "./data/f_and_sx_bnew_%.1f_g.txt" b
	b1 = load(path1)
end

# ╔═╡ 2cb83a03-03b2-4958-bcc0-4877b6ad658d
gamma = [i for i in range(0,5,step = 0.2)]

# ╔═╡ 9d62aefc-6930-46d6-8a68-225c25be3d51
z = pauli('z')

# ╔═╡ c7d1df2b-4261-4347-a619-778a5f13ecbf
num = 10

# ╔═╡ b98faf74-3372-4f36-b9e9-2318355ecaa9
begin
	ge2 = zeros(length(beta), num); ge3 = zeros(length(beta), num)
	dge2 = zeros(length(beta), num-1); dge3 = zeros(length(beta), num-1)
	dge2_beta = zeros(length(beta), num-1); dge3_beta = zeros(length(beta), num-1)
	for i = 1:length(beta)
		β = beta[i]; key = string(beta[i])
		ψ = cmps(g1[key][2][:,:,1], g1[key][2][:,:,1])
		gk2 = ψ * ψ; ge2[i,:] = eigvals(gk2)[1:num]
		dge2[i,:] = ge2[i,2:end] .- ge2[i,1]
		dge2_beta[i,:] = dge2[i,:] * beta[i]
		gk3 = ψ * w1 * ψ; ge3[i,:] = eigvals(gk3)[1:num]
		dge3[i,:] = ge3[i,2:end] .- ge3[i,1]
		dge3_beta[i,:] = dge3[i,:] * beta[i]
	end	
end

# ╔═╡ 16d87d8a-d159-4ec6-8969-1acdaf787341
begin
	# plot gk2
	gti = @sprintf "Γ/J = %.1f, eigen values of |--|" g
	plot(beta, dge2_beta, ls=:dash, shape=:circle)
	plot!(title= gti, xlabel = "β", ylabel = "ΔE/T",legend=:best)
end

# ╔═╡ 9f64825a-6432-4987-b4c8-abf4990664c7
begin
	# plot gk2
	gti2 = @sprintf "Γ/J = %.1f, eigen values of |-+-|" g
	x = log.(T)
	plot(x, dge3_beta, ls=:dash, shape=:circle)
	#plot!(x, -0.75 .* x,lw = 1.5, c=:black, label=("-0.75log(T)"))
	plot!(title= gti2, xlabel = "log(T)", ylabel = "ΔE/T",legend=:topright)
end

# ╔═╡ 4572629a-e1dd-4d48-9a62-6e130b25cea9
begin
	e2 = zeros(length(gamma),num); e3 = zeros(length(gamma),num)
	de2 = zeros(length(gamma), num-1); de3 = zeros(length(gamma), num-1)
	for i = 1:length(gamma)
		gi = gamma[i]; key = string(gamma[i])
		w = TFIsing(1.0, gi)
		ψi = cmps(b1[key][2][:,:,1], b1[key][2][:,:,1])
		k2 = ψi * ψi; e2[i,:] = eigvals(k2)[1:num]
		de2[i,:] = e2[i,2:end] .- e2[i,1]
		k3 = ψi * w * ψi; e3[i,:] = eigvals(k3)[1:num]
		de3[i,:] = e3[i,2:end] .- e3[i,1]
	end	
end

# ╔═╡ b236eed5-1ed8-4fee-9d26-0561ad067935
begin
	# plot k2
	bti = @sprintf "β=%.1f, eigen values of |--|" b
	plot(gamma, de2, ls=:dash, shape=:circle)
	plot!(title=bti, xlabel = "Γ/J", ylabel = "ΔE = Ei-E1",legend=:topleft)
end

# ╔═╡ 27300b00-fb08-43c3-b691-e3d65973e18b
begin
	# plot k3
	bti2 = @sprintf "β=%.1f, eigen values of |-+-|" b
	plot(gamma, de3, ls=:dash, shape=:circle)
	plot!(gamma, 4 .* abs.(gamma .- 1.0),lw = 1.5, label=("2|Δ|"))
	plot!(gamma, 8 .* abs.(gamma .- 1.0),lw = 1.5, label=("4|Δ|"))
	plot!(title= bti2, xlabel = "Γ/J", ylabel = "ΔE = Ei-E1",legend=:topleft)
end

# ╔═╡ 0065491c-a633-11eb-25e8-291ac2c6fa6e
pwd()

# ╔═╡ d8d0e971-fb53-47dd-8dc5-5808e54bc8af
begin
	#compare free energy
	f1 = readdlm(p1)
	f2 = readdlm(p2)
	f3 = readdlm(p3)
end

# ╔═╡ 443aa1dc-e4fa-4766-a51d-ea4ddcba3494
df1 = f1[:,3] - f1[:,2]

# ╔═╡ 80c3bc9e-d417-4682-861d-db0fef01b86c
df2 = f2[:,3] - f2[:,2]

# ╔═╡ ef56a38e-15a8-4aee-a870-c11171711e26
df3 = f3[:,3] - f3[:,2]

# ╔═╡ 10b10dc3-4cae-4b2b-a725-a760abdbe7e1
begin
	tf=@sprintf "free energy loss compare of β=%.1f with different ψ0" b
	plot(gamma, df1, m=:circle, lw = 2, label="+-| init")
	plot!(gamma, df2, m=:circle, lw = 2, label="random init")
	plot!(gamma, df3, m=:circle, lw = 2, label="previous init")
	plot!(title=tf, xlabel="Γ/J", ylabel="free energy loss",legend=:best)
end

# ╔═╡ 58a5077f-0cab-48e4-b2a3-c0e7375615f6
begin
	plot(gamma, df1, m=:circle, lw = 2, label="+-| init")
	plot!(gamma, df2, m=:circle, lw = 2, label="random init")
	plot!(title=tf, xlabel="Γ/J", ylabel="free energy loss", ylims=(-1.e-10,6.e-9), legend=:best)
end

# ╔═╡ Cell order:
# ╠═a63b4ad6-46e6-4635-8037-c24b20abf070
# ╠═181610df-4507-4cec-a7e0-2f919aec32bd
# ╟─da8e7582-314f-401b-9406-a334c299c6b2
# ╟─ddaf571f-8bea-4c47-b8f7-1d7499de4719
# ╟─0cd9aa6f-c170-4fd1-9284-e431b76cd821
# ╠═b98faf74-3372-4f36-b9e9-2318355ecaa9
# ╠═16d87d8a-d159-4ec6-8969-1acdaf787341
# ╠═9f64825a-6432-4987-b4c8-abf4990664c7
# ╠═b6015a68-602b-4883-9cc3-9e57c8dffdd7
# ╠═2cb83a03-03b2-4958-bcc0-4877b6ad658d
# ╠═4572629a-e1dd-4d48-9a62-6e130b25cea9
# ╠═b236eed5-1ed8-4fee-9d26-0561ad067935
# ╠═27300b00-fb08-43c3-b691-e3d65973e18b
# ╠═9d62aefc-6930-46d6-8a68-225c25be3d51
# ╠═c7d1df2b-4261-4347-a619-778a5f13ecbf
# ╠═91909ae8-cb1f-4653-9ba6-144bd4d01747
# ╠═76d940bc-09c4-422c-ab5a-50c8be7534e2
# ╠═86574854-9b9c-4fa8-940f-20446705fce6
# ╠═0065491c-a633-11eb-25e8-291ac2c6fa6e
# ╠═d8d0e971-fb53-47dd-8dc5-5808e54bc8af
# ╠═443aa1dc-e4fa-4766-a51d-ea4ddcba3494
# ╠═80c3bc9e-d417-4682-861d-db0fef01b86c
# ╠═ef56a38e-15a8-4aee-a870-c11171711e26
# ╠═10b10dc3-4cae-4b2b-a725-a760abdbe7e1
# ╠═58a5077f-0cab-48e4-b2a3-c0e7375615f6
