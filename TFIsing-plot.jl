### A Pluto.jl notebook ###
# v0.14.7

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

# ╔═╡ bac5d06e-17d2-44a3-b7b5-09f924592673
begin
	# gamma = [0.1, 0.5, 0.9, 1, 1.1, 2, 4]
	g_2 = 1.0
	path3 = @sprintf "./data/gnew_%.1f.jld" g_2
	g2 = load(path3)
end

# ╔═╡ a7c05b5c-d77d-4d71-9a5f-139a621968d8
w2 = TFIsing(1.0, g_2)

# ╔═╡ b6d66d89-6a5d-47ab-bd89-fa8e0e54fbf0
begin
	# Fit the lowest energy level (E2-E1)
	model(t,p) = p[1] * t .+ p[2]
	p0 = [-0.75, 0]
end

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
gamma = [i for i in range(0,2,step = 0.2)]

# ╔═╡ aa944832-e966-4d56-9a53-42a92b8dd75b
# v = [gamma, 128, 128] eigen vectors, length(gamma)=26

# ╔═╡ 9d62aefc-6930-46d6-8a68-225c25be3d51
z = pauli('z')

# ╔═╡ ed2f430a-5a9d-47ac-ac3a-fa07ff6a7680
begin
	eye = Matrix(1.0I, 8, 8)
	sz = eye ⊗ z ⊗ eye
end

# ╔═╡ c7d1df2b-4261-4347-a619-778a5f13ecbf
num = 10

# ╔═╡ b98faf74-3372-4f36-b9e9-2318355ecaa9
begin
	ge2 = zeros(length(beta), num); ge3 = zeros(length(beta), num)
	dge2 = zeros(length(beta), num-1); dge3 = zeros(length(beta), num-1)
	dge2_beta = zeros(length(beta), num-1); dge3_beta = zeros(length(beta), num-1)
	for i = 1:length(beta)
		β = beta[i]; key = string(beta[i])
		ψ = cmps(g1[key][2][:,:,1], g1[key][2][:,:,2])
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
	plot!(x, -0.75 .* x,lw = 1.5, c=:black, label=("-0.75log(T)"))
	plot!(title= gti2, xlabel = "log(T)", ylabel = "ΔE/T",legend=:topright)
end

# ╔═╡ f01991c3-5906-440f-9153-a4d3fa64e3c0
begin
	ge_2 = zeros(length(beta), num); ge_3 = zeros(length(beta), num)
	dge_2 = zeros(length(beta), num-1); dge_3 = zeros(length(beta), num-1)
	dge_2_beta = zeros(length(beta), num-1); dge_3_beta = zeros(length(beta), num-1)
	v2 = zeros(length(beta), 128, 128)
	for i = 1:length(beta)
		β = beta[i]; key = string(beta[i])
		ψ = tocmps(g2[key][2])
		gk2 = ψ * ψ; ge_2[i,:] = eigvals(gk2)[1:num]
		dge_2[i,:] = ge2[i,2:end] .- ge2[i,1]
		dge_2_beta[i,:] = dge2[i,:] * beta[i]
		gk3 = ψ * w2 * ψ
		e, v2[i,:,:] = eigen(gk3)
		ge_3[i,:] = e[1:num]
		dge_3[i,:] = ge_3[i,2:end] .- ge_3[i,1]
		dge_3_beta[i,:] = dge_3[i,:] * beta[i]
	end	
end

# ╔═╡ 6641f603-2764-49cd-8d20-c3aeabc09446
begin
	sz2 = zeros(length(beta),128,128)
	for i = 1:length(beta)
		sz2[i,:,:] = v2[i,:,:]' * sz * v2[i,:,:]
	end
end

# ╔═╡ 84514442-d9d7-49bd-8394-94b7c551976d
begin
	# sz[i,i] * sz[i,i] delta(0) part
	# gamma = 1: weak linear 1/T dependence
	y3 = zeros(length(beta), num)
	for i = 1: length(beta), j = 1:num
		y3[i,j] = sz2[i,j,j] * sz2[i,j,j]
	end
	# plot 
	yti3 = @sprintf "Γ/J =%.1f, sz[i,i] * sz[i,i]" g_2
	plot(x, y3, ls=:dash, shape=:circle)
	plot!(title= yti3, xlabel = "log(T)", ylabel = "|z_ii|^2 ", legend=:bottomleft)
end

# ╔═╡ 831f2394-3fd1-45b7-9cbc-32886bb10b93
begin # |sz[i,i]|^2 * e^{-Em/T}
	yd3_beta = zeros(length(beta), num-1)
	for i = 1: length(beta), j = 1:num-1
		yd3_beta[i,j] = -log(y3[i,j+1]) + dge3_beta[i,j]
	end
	plot(x, yd3_beta, ls=:dash, shape=:circle)
	plot!(x, -0.75 .* x,lw = 1.5, c=:black, label=("-0.75log(T)"))
	plot!(xlabel = "log(T)", ylabel = "log(|z_ii|^2)+ΔE/T ", legend=:best)
end

# ╔═╡ e51012e1-832e-4c6d-b4b3-db36be55f7c0
begin
	x1 = x[10:end] # inear part
	ydata1 = dge3_beta[10:end,1]
	ydata2 = yd3_beta[10:end,1]
	fit1 = curve_fit(model, x1, ydata1,p0); param1 = fit1.param
	lab1 = @sprintf "fit ΔE1/T = %.2f * log(T) + %.2f" param1[1] param1[2]
	fit2 = curve_fit(model, x1, ydata2,p0); param2 = fit2.param
	lab2 = @sprintf "fit log|Szii|^2 e^{-ΔE1/T} = %.2f * log(T) + %.2f" param1[1] param1[2]
	nydata1 = model(x1,param1)
	nydata2 = model(x1,param2)
	plot(x1, ydata1, ls=:dash, shape=:circle, label="ΔE1/T")
	plot!(x1, nydata1, lw=2, label="fit ΔE1/T")
end

# ╔═╡ ddb84c1e-15f0-4d42-b756-f63da6964866
fit2.param

# ╔═╡ 2a448062-f5d2-4422-b8f0-b3c070c023e1
begin
	# sz[i,j] * sz[j,i]
	y4 = zeros(length(beta), num)
	for i = 1: length(beta), j = 1:num
		y4[i,j] = sz2[i,j,:]' * sz2[i,:,j] - 1
	end
	# plot 
	yti4 = @sprintf "Γ/J =%.1f, sz[i,j] * sz[j,i]-1" g_2
	plot(beta, y4, ls=:dash, shape=:circle)
	plot!(title= yti4, xlabel = "β",legend=:best)
end

# ╔═╡ 4572629a-e1dd-4d48-9a62-6e130b25cea9
begin
	e2 = zeros(length(gamma),num); e3 = zeros(length(gamma),num)
	de2 = zeros(length(gamma), num-1); de3 = zeros(length(gamma), num-1)
	v = zeros(length(gamma),128, 128)
	for i = 1:length(gamma)
		gi = gamma[i]; key = string(gamma[i])
		w = TFIsing(1.0, gi)
		ψi = cmps(b1[key][2][:,:,1], b1[key][2][:,:,2])
		k2 = ψi * ψi; e2[i,:] = eigvals(k2)[1:num]
		de2[i,:] = e2[i,2:end] .- e2[i,1]
		k3 = ψi * w * ψi
		e, v[i,:,:] = eigen(k3)
		e3[i,:] = e[1:num] # first num eigen values
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
	plot!(gamma, 2 .* abs.(gamma .- 1.0),lw = 1.5, label=("|Δ|"))
	plot!(gamma, 4 .* abs.(gamma .- 1.0),lw = 1.5, label=("2|Δ|"))
	plot!(title= bti2, xlabel = "Γ/J", ylabel = "ΔE = Ei-E1",legend=:topleft)
end

# ╔═╡ 61cd867d-b295-4d1f-b00c-df224c7fc1de
begin
	sz1 = zeros(length(gamma),128,128)
	for i = 1:length(gamma)
		sz1[i,:,:] = v[i,:,:]' * sz * v[i,:,:]
	end
end

# ╔═╡ 8bb6c29a-7e45-492d-9aaf-ead80b16e854
begin
	# sz[i,i] * sz[i,i] delta(0) part
	y1 = zeros(length(gamma), num)
	for i = 1: length(gamma), j = 1:num
		y1[i,j] = sz1[i,j,j] * sz1[i,j,j]
	end
	# plot 
	yti1 = @sprintf "β=%.1f, sz[i,i] * sz[i,i] delta(0) part" b
	plot(gamma, y1, ls=:dash, shape=:circle)
	plot!(title= yti1, xlabel = "Γ/J", legend=:bottomleft)
end

# ╔═╡ fa127125-f8a7-4af1-90cd-6d21ba754db7
begin
	# sz[i,i+1] * sz[i+1,i] delta(0) part
	y2 = zeros(length(gamma), num)
	for i = 1: length(gamma), j = 1:num
		y2[i,j] += sz1[i,j,:]' * sz1[i,:,j] - 1
	end
	# plot 
	yti2 = @sprintf "β=%.1f, sz[i,j] * sz[j,i] - 1" b
	plot(gamma, y2, ls=:dash, shape=:circle)
	plot!(title= yti2, xlabel = "Γ/J",legend=:bottomleft)
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
# ╠═181610df-4507-4cec-a7e0-2f919aec32bd
# ╠═da8e7582-314f-401b-9406-a334c299c6b2
# ╟─ddaf571f-8bea-4c47-b8f7-1d7499de4719
# ╟─0cd9aa6f-c170-4fd1-9284-e431b76cd821
# ╠═b98faf74-3372-4f36-b9e9-2318355ecaa9
# ╠═16d87d8a-d159-4ec6-8969-1acdaf787341
# ╠═9f64825a-6432-4987-b4c8-abf4990664c7
# ╠═bac5d06e-17d2-44a3-b7b5-09f924592673
# ╟─a7c05b5c-d77d-4d71-9a5f-139a621968d8
# ╠═f01991c3-5906-440f-9153-a4d3fa64e3c0
# ╠═6641f603-2764-49cd-8d20-c3aeabc09446
# ╠═84514442-d9d7-49bd-8394-94b7c551976d
# ╠═831f2394-3fd1-45b7-9cbc-32886bb10b93
# ╠═b6d66d89-6a5d-47ab-bd89-fa8e0e54fbf0
# ╠═e51012e1-832e-4c6d-b4b3-db36be55f7c0
# ╠═ddb84c1e-15f0-4d42-b756-f63da6964866
# ╠═2a448062-f5d2-4422-b8f0-b3c070c023e1
# ╠═b6015a68-602b-4883-9cc3-9e57c8dffdd7
# ╠═2cb83a03-03b2-4958-bcc0-4877b6ad658d
# ╠═4572629a-e1dd-4d48-9a62-6e130b25cea9
# ╠═b236eed5-1ed8-4fee-9d26-0561ad067935
# ╠═27300b00-fb08-43c3-b691-e3d65973e18b
# ╠═aa944832-e966-4d56-9a53-42a92b8dd75b
# ╠═61cd867d-b295-4d1f-b00c-df224c7fc1de
# ╠═8bb6c29a-7e45-492d-9aaf-ead80b16e854
# ╠═fa127125-f8a7-4af1-90cd-6d21ba754db7
# ╠═9d62aefc-6930-46d6-8a68-225c25be3d51
# ╠═ed2f430a-5a9d-47ac-ac3a-fa07ff6a7680
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
