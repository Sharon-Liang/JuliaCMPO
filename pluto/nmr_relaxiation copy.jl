### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 9f3511fc-e24d-11eb-130d-15fd034ec045
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

# ╔═╡ 838f933c-db31-4887-b206-562bd060746f
#gamma = [0.1, 0.5, 1.0, 2.0, 4.0, 8.0]
g = 4.0

# ╔═╡ f4fda113-17ce-4312-93ff-fb361815b970
begin
	path = @sprintf "../data/g_%.1f.jld" g
	data = load(path)
	"load data"
end

# ╔═╡ c3f0d275-35f2-4c35-a8b4-ffab0cb07732
md"""
β = [i for i in range(1,20,step=0.02)]
"""

# ╔═╡ 3f990c4b-97da-4336-8e34-1bd496162328
lowT = [i for i in range(10, 20, step=0.04)]

# ╔═╡ c30979b4-97dc-4252-ad32-eb8ead672c8a
highT = [i for i in range(1, 20, step=0.2)] 

# ╔═╡ 9dd8334d-7664-44f9-b052-b3858d913066
function linearfit(xdata::Vector, ydata::Vector)
	p0 = zeros(2)
	if g==1.
		@. model(x,p) = p[1]*x + p[2] 
		fit= curve_fit(model, xdata, ydata, p0)
		fname = @sprintf "%.3f x + %.3f" fit.param[1] fit.param[2]
		f = x -> model(x, fit.param)
	elseif g != 1.
		Δ=abs(1-g)*2
		@. model1(x,p) = p[1]*Δ*x + p[2] 
		fit= curve_fit(model1, xdata, ydata, p0)
		fname = @sprintf "%.3f|Δ|x + %.3f" fit.param[1] fit.param[2]
		f = x -> model(x, fit.param)
	end
	return f, fname
end


# ╔═╡ 31b2e157-52ac-4b56-a3bd-fc0c18c97499
begin
	function expfunc(p::Number)
		@. f(x) = exp(p*x)
		fname = @sprintf "exp(%.1f x)" p
		return f, fname
	end
	
	x1 = [i for i in range(1,20, length=100)]
	f1, fname1 = expfunc(1)
	f2, fname2 = expfunc(-1)
end

# ╔═╡ db9bfd32-557a-4e77-82ee-c418ef07fa22
begin
	plot(x1, f1(x1), lw= 2, label = fname1)
	#plot!(yaxis = :log)
end

# ╔═╡ 66f97cd2-0488-4721-bc64-d9277421292f
begin
	plot(1 ./x1, f1(x1), lw= 2, label = fname1)
	plot!(xlabel = "1/x")
	#plot!(yaxis = :log)
end

# ╔═╡ c661e404-a709-4cf2-a344-f02d4d4077be
begin
	plot(x1, f2(x1), lw = 2, label = fname2)
	#plot!(yaxis = :log)
end

# ╔═╡ 68ac54e2-79ad-4a0b-b0ef-5ff8b3f28e39
begin
	plot(1 ./ x1, f2(x1), lw = 2, label = fname2)
	plot!(xlabel = "1/x")
end

# ╔═╡ c3b1983c-f32a-401c-9720-c53bbd0a7825
begin
	struct nmr_result
		fig
		xdata::Vector{Float64}
		ydata::Vector{Float64}
		linefit::Bool
	end
	"nmr_result struct"
end

# ╔═╡ 05f4d459-5bc7-4dbf-ab13-42041c0cf8e7
begin
	z = make_operator(pauli(:z),8)
	x = make_operator(pauli(:x),8)
	"x, z operator"
end

# ╔═╡ dfe621e5-a860-47e0-b670-2fd4f57eae68
function nmr_relaxiation(beta::Vector,η::Float64)
	w = TFIsing(1., g)
	len = length(beta)
	res = zeros(len)
	for i = 1:len
		β = beta[i]; key = string(β)
		ψ = tocmps(data[key][2])
		res[i] = structure_factor(0,z,z,ψ,w,β,η=η,method=:S)
	end
	return res
end

# ╔═╡ 4dac8cf2-8e9d-4aa5-9061-1228332a08b9
function nmr_relaxiation(beta::Vector)
	w = TFIsing(1., g)
	len = length(beta)
	res = zeros(len)
	for i = 1:len
		β = beta[i]; key = string(β)
		ψ = tocmps(data[key][2])
		res[i] = π/4*β*correlation_2time(β/2,z,z,ψ,w,β)
	end
	return res
end

# ╔═╡ 5f036037-0e4a-42a8-be10-2e467ca1015d
function plot_nmr(beta::Vector, η::Float64; 
		fit::Bool=false, xtick::Symbol=:beta,
		xscale::Symbol=:plain, yscale::Symbol=:plain )
	ptitle = @sprintf "Γ/J=%.1f" g
	datalab = @sprintf "η=%.3f" η
	ydata = nmr_relaxiation(beta, η)

	if xtick == :beta
		xdata = beta
		xscale == :log ? xlab = "log(β)" : xlab = "β"
	elseif xtick == :T
		xdata = 1 ./beta
		xscale == :log ? xlab = "log(T)" : xlab = "T"
	else
		@error "xtick should be :beta or :T"
	end
	
	xscale == :log ? xdata = log.(xdata) : xdata = xdata
	
	yscale == :log ? ydata = log.(ydata) : ydata = ydata
	
	yscale == :log ? ylab = "log(1/T1)" : ylab = "1/T1"
	
	fig = plot(xdata, ydata, line=(:dash), marker=(:circle,3),label=datalab)
	
	if g == 1
		linefit = (yscale == :log) && (xscale == :log)
	else
		linefit = (yscale == :log) && (xtick == :beta) && (xscale == :plain)
	end
	
	plot!(title = ptitle, xlabel=xlab, ylabel = ylab)
	return nmr_result(fig, xdata, ydata, linefit)
end

# ╔═╡ e4f3cb7b-b191-4e0d-95ff-c1c5ccd51f1d
function plot_nmr(beta::Vector; 
		fit::Bool=false, xtick::Symbol=:beta,
		xscale::Symbol=:plain, yscale::Symbol=:plain )
	ptitle = @sprintf "Γ/J=%.1f" g
	datalab = @sprintf "β/2 average" 
	ydata = nmr_relaxiation(beta)

	if xtick == :beta
		xdata = beta
		xscale == :log ? xlab = "log(β)" : xlab = "β"
	elseif xtick == :T
		xdata = 1 ./beta
		xscale == :log ? xlab = "log(T)" : xlab = "T"
	else
		@error "xtick should be :beta or :T"
	end
	
	xscale == :log ? xdata = log.(xdata) : xdata = xdata
	
	yscale == :log ? ydata = log.(ydata) : ydata = ydata
	
	yscale == :log ? ylab = "log(1/T1)" : ylab = "1/T1"
	
	fig = plot(xdata, ydata, line=(:dash), marker=(:circle,3),label=datalab)
	
	if g == 1
		linefit = (yscale == :log) && (xscale == :log)
	else
		linefit = (yscale == :log) && (xtick == :beta) && (xscale == :plain)
	end
	
	plot!(title = ptitle, xlabel=xlab, ylabel = ylab)
	return nmr_result(fig, xdata, ydata, linefit)
end

# ╔═╡ 6028341a-c20b-46bb-b6d9-d654eada7d10
begin
	p1= plot_nmr(highT, 0.05, xtick = :beta, yscale=:log)
	p1.fig
	#plot!(ylim=(minimum(p1.ydata)*1.1, maximum(p1.ydata)*0.9 ))
end

# ╔═╡ b07be761-27f7-4b9e-83ce-d0b7891cfd85
p1.linefit

# ╔═╡ 5134ad41-e057-4fff-a41f-d5f3a290c3f3
begin
	p2= plot_nmr(lowT, xtick = :beta, yscale=:log)
	p2.fig
end

# ╔═╡ 6a4197dd-e1a0-4496-b389-790b24692b0f
p2.linefit

# ╔═╡ 6b8708f3-6504-423b-b962-ceb69f0551a3
pwd()

# ╔═╡ Cell order:
# ╟─838f933c-db31-4887-b206-562bd060746f
# ╟─f4fda113-17ce-4312-93ff-fb361815b970
# ╟─c3f0d275-35f2-4c35-a8b4-ffab0cb07732
# ╟─3f990c4b-97da-4336-8e34-1bd496162328
# ╟─c30979b4-97dc-4252-ad32-eb8ead672c8a
# ╟─dfe621e5-a860-47e0-b670-2fd4f57eae68
# ╟─4dac8cf2-8e9d-4aa5-9061-1228332a08b9
# ╟─9dd8334d-7664-44f9-b052-b3858d913066
# ╟─5f036037-0e4a-42a8-be10-2e467ca1015d
# ╟─e4f3cb7b-b191-4e0d-95ff-c1c5ccd51f1d
# ╠═6028341a-c20b-46bb-b6d9-d654eada7d10
# ╠═b07be761-27f7-4b9e-83ce-d0b7891cfd85
# ╠═5134ad41-e057-4fff-a41f-d5f3a290c3f3
# ╠═6a4197dd-e1a0-4496-b389-790b24692b0f
# ╟─31b2e157-52ac-4b56-a3bd-fc0c18c97499
# ╟─db9bfd32-557a-4e77-82ee-c418ef07fa22
# ╟─66f97cd2-0488-4721-bc64-d9277421292f
# ╟─c661e404-a709-4cf2-a344-f02d4d4077be
# ╟─68ac54e2-79ad-4a0b-b0ef-5ff8b3f28e39
# ╠═c3b1983c-f32a-401c-9720-c53bbd0a7825
# ╟─05f4d459-5bc7-4dbf-ab13-42041c0cf8e7
# ╟─9f3511fc-e24d-11eb-130d-15fd034ec045
# ╟─6b8708f3-6504-423b-b962-ceb69f0551a3
