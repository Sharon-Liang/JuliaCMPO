### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ f3e91cc8-41cc-11ec-0cc4-0fe7d3ecfd69
begin
	using OMEinsum
	using Zygote, FiniteDiff, ForwardDiff
	using LinearAlgebra, GenericLinearAlgebra
	using StatsFuns
	"load packages"
end

# ╔═╡ 56a94c71-1329-40cd-8bed-3748307ba2b4
md"""
### free energy
"""

# ╔═╡ 22066571-3e2b-43a1-862c-6a40bc3f4b47
md"""
### product
"""

# ╔═╡ 336ee5a5-eee2-416d-a72e-27cdff89369c
md"""
### eigval
"""

# ╔═╡ 36dea9c1-f358-4bfc-b243-2c33919ae5c6
md"""
### logsumexp
"""

# ╔═╡ 9bc2dd49-ee54-4e58-bac0-13426d9a7098
md"""
### logtrexp
"""

# ╔═╡ 769b2024-6459-44c2-ab36-01b74944b921
begin
	struct CMPS
		Q::Array{<:Number}
		R::Array{<:Number}
	end
		
	struct CMPO
		Q::Array{<:Number}  # onsite
		R::Array{<:Number}  # interaction, column vector
		L::Array{<:Number}  # interaction, row vector
		P::Array{<:Number}  # long-range
	end
	"define struct"
end

# ╔═╡ 88ae1058-abcf-4567-86dd-8e1f1c5f7871
begin	
	function toarray(ψ::CMPS)
	    sq = size(ψ.Q)
	    sr = size(ψ.R) #sr: dimension of ψ.R array
	    if length(sr) == 2 #sr = 2, ψ.R is a matrix
	        Q = reshape(ψ.Q, sq[1],sq[2],1)
	        R = reshape(ψ.R, sr[1],sr[2],1)
	    elseif length(sr) > 2
	        Q = reshape(ψ.Q, sq[1],sq[2],1)
	        R = ψ.R
	    else
	        @error "Improper CMPS"
	    end
	    return cat(Q,R,dims=3)
	end
	
	function tocmps(A::Array{<:Number,3})
	    d = size(A)[3]
	    if d == 2
	        return CMPS(A[:,:,1],A[:,:,2])
	    else
	        return CMPS(A[:,:,1],A[:,:,2:end])
	    end
	end
	"conversion of cmps and array"
end

# ╔═╡ 2bbca5ed-12bc-4646-a401-d691e9c138eb
begin
	import Base: *
	function ⊗(A::Matrix{<:Number}, B::Matrix{<:Number})
	    (r1, c1) = size(A)
	    (r2,c2) = size(B)
	    return reshape(ein"ij,kl->kilj"(A, B), r1*r2, c1*c2)
	end
	
	function ⊗(A::Matrix{<:Number}, B::Array{<:Number,3})
	    (r1,c1) = size(A)
	    (r2,c2,d) = size(B)
	    return reshape(ein"ij,klm->kiljm"(A, B), r1*r2, c1*c2, d)
	end
	
	function ⊗(A::Array{<:Number,3}, B::Matrix{<:Number})
	    (r1,c1,d) = size(A)
	    (r2,c2) = size(B)
	    return reshape(ein"ijm,kl->kiljm"(A, B), r1*r2, c1*c2, d)
	end
	
	function ⊗(A::Array{<:Number,3}, B::Array{<:Number,3})
	    (r1,c1,d1) = size(A)
	    (r2,c2,d2) = size(B)
	    if d1 != d2 @error "Dimension mismatch!" end
	    return reshape(ein"ijm,klm->kilj"(A, B), r1*r2, c1*c2)
	end
	
	function ⊗(A::Array{<:Number,4}, B::Array{<:Number,3})
	    (r1,c1,d1,f) = size(A)
	    (r2,c2,d2) = size(B)
	    if f != d2 @error "Dimension mismatch!" end
	    return reshape(ein"ijnm,klm->kiljn"(A, B), r1*r2, c1*c2, d1)    
	end
	
	function ⊗(A::Array{<:Number,3}, B::Array{<:Number,4}) 
	    (r1,c1,d1) = size(A)
	    (r2,c2,d2,f) = size(B)
	    if d1 != d2 @error "Dimension mismatch!" end
	    return reshape(ein"ijm,klmn->kiljn"(A, B), r1*r2, c1*c2, f)    
	end
	
	function ⊗(A::Array{<:Number,4}, B::Array{<:Number,4})
	    (r1,c1,d1,f1) = size(A)
	    (r2,c2,d2,f2) = size(B)
	    if f1 != d2 @error "Dimension mismatch!" end
	    return reshape(ein"ijpm,klmq->kiljpq"(A, B), r1*r2, c1*c2, d1, f2)    
	end
	
	function *(sl::CMPS, sr::CMPS)
	    li = Matrix(1.0I,size(sl.Q))
	    ri = Matrix(1.0I,size(sr.Q))
	    K = li ⊗ sr.Q + sl.Q ⊗ ri + sl.R ⊗ sr.R
	    return -K
	end
	
	function *(o::CMPO, s::CMPS)
	    oi = Matrix(1.0I,size(o.Q))
	    si = Matrix(1.0I,size(s.Q))
	    Q = oi ⊗ s.Q + o.Q ⊗ si + o.L ⊗ s.R
	    R = o.R ⊗ si + o.P ⊗ s.R
	    return CMPS(Q, R)
	end
	
	function *(s::CMPS, o::CMPO)
	    oi = Matrix(1.0I,size(o.Q))
	    si = Matrix(1.0I,size(s.Q))
	    Q = s.Q ⊗ oi + si ⊗ o.Q + s.R ⊗ o.R
	    R = si ⊗ o.L + s.R ⊗ o.P
	    return CMPS(Q, R)
	end
	"define products of cmps and cmpo"
end

# ╔═╡ 4de5fbf3-afce-4d12-883d-ab3b3b6da1b6
begin
	function symmetrize(A::T where T<:AbstractMatrix)
	    (A + A')/2
	end
	
	function logtrexp(A::T where T<:AbstractMatrix)
	    #if ishermitian(A) == false
	    if isapprox(A,A') == false
	        error("The input matrix should be hermitian")
	    end
	    A = symmetrize(A) |> Hermitian
	    return eigvals(A) |> logsumexp
	end
	
	function free_energy(ψ::CMPS, W::CMPO, β::Real)
	    K = ψ * W * ψ |> symmetrize |> Hermitian
	    H = ψ * ψ |> symmetrize |> Hermitian
	    res = logtrexp(-β*K)- logtrexp(-β*H)
	    return -1/β * res
	end
	
	function free_energy(param::Array{<:Number,3}, W::CMPO, β::Real)
	    free_energy(tocmps(param), W, β)
	end
	"free energy function"
end

# ╔═╡ d5b61aa8-5125-4a3d-804f-26b55879c0f5
begin
	function K(ψ::CMPS, W::CMPO, β::Real)
	    #K = ψ * W * ψ |> symmetrize |> Hermitian
		K = ψ * ψ |> symmetrize |> Hermitian
	    return sum(K)
	end
	
	function K(param::Array{<:Number,3}, W::CMPO, β::Real)
	    K(tocmps(param), W, β)
	end
end

# ╔═╡ ad7aa9c3-d546-4e1e-931a-d56377532cf2
begin
	function Keig(ψ::CMPS, W::CMPO, β::Real)
	    #K = ψ * W * ψ |> symmetrize |> Hermitian
		K = ψ * ψ |> symmetrize |> Hermitian
	    return eigvals(K) |> sum
	end
	
	function Keig(param::Array{<:Number,3}, W::CMPO, β::Real)
	    Keig(tocmps(param), W, β)
	end
end

# ╔═╡ 123622a5-f727-4426-9c6e-5fdf749ec610
begin
	function Keiglog(ψ::CMPS, W::CMPO, β::Real)
		#K = ψ * W * ψ |> symmetrize |> Hermitian
		K = ψ * ψ |> symmetrize |> Hermitian
		return eigvals(K) |> logsumexp
	end
		
	function Keiglog(param::Array{<:Number,3}, W::CMPO, β::Real)
		Keig(tocmps(param), W, β)
	end
end

# ╔═╡ f50f2e48-3da1-49e5-9f2e-c3eb22e3f6b4
begin
	function logtr(ψ::CMPS, W::CMPO, β::Real)
		#K = ψ * W * ψ |> symmetrize |> Hermitian
		K = ψ * ψ |> symmetrize |> Hermitian
		return logtrexp(-β*K)
	end
		
	function logtr(param::Array{<:Number,3}, W::CMPO, β::Real)
		logtr(tocmps(param), W, β)
	end
end

# ╔═╡ 91fc3a50-0420-4e3f-8e15-58fe446ce0be
begin
	function pauli(symbol::Symbol)
	    if symbol==:x return [0. 1.; 1. 0.]
	    elseif symbol==:y return [0. -1im; 1im 0.]
	    elseif symbol==:z return [1. 0.; 0. -1.]
	    elseif symbol==:+ return [0. 1.; 0. 0.]
	    elseif symbol==:- return [0. 0.; 1. 0.]
	    else
	        error("The input should be :x,:y,:z,:+,:-.")
	    end
	end
	
	function TFIsing(J::Real, Γ::Real; field::Symbol=:N, η::Float64 = 1.e-2)
	    if field == :N
	        h = zeros(2,2)
	    else
	        h = η .* pauli(field)
	    end
	    return CMPO(Γ*pauli(:x)+h, √J*pauli(:z), √J*pauli(:z), zeros(2,2))
	end
	
	function init_cmps(χ::Int64, W::CMPO)
	    # r = 0 case
	    d = size(W.Q)[1];  (q,r) = divrem(log(d,χ), 1)
	    ψ = CMPS(W.Q, W.R)
	    if r == 0
	        for i = 1:Integer(q-1)  ψ = W * ψ  end
	    else
	        error("Not support yet :)")
	    end
	    return ψ
	end
	"utilities"
end

# ╔═╡ fc634506-632f-47ba-8a07-b5ce108c2012
w = TFIsing(1.0,1.0)

# ╔═╡ fe3764e3-f309-41f5-8d2d-8d1783d87cba
test1 = x->free_energy(reshape(x,2,2,2),w,20)

# ╔═╡ 97e58e64-7dfc-4414-a7e7-8d1f38e9a4ec
test2 = x->K(reshape(x,2,2,2),w,20)

# ╔═╡ 10516ffa-6e74-4d95-9ebd-492166e8ca28
test3 = x->Keig(reshape(x,2,2,2),w,20)

# ╔═╡ 0208996b-a133-4c72-9a5f-048f963924ac
test4 = x->Keiglog(reshape(x,2,2,2),w,20)

# ╔═╡ 23065c11-ebe1-47a2-a770-c4283dd2c36b
test5 = x->logtr(reshape(x,2,2,2),w,20)

# ╔═╡ 40e6e0a9-ca69-47ba-b982-ea2e9d098578
v = init_cmps(2,w) |> toarray |> vec

# ╔═╡ 604dfba0-5d9d-4fbf-82bc-7f19f4032c5b
begin
	g1 = zeros(length(v),3)
	g1[:,1] = ForwardDiff.gradient(test1, v)
	g1[:,2] = Zygote.gradient(test1, v)[1]
	g1[:,3] = FiniteDiff.finite_difference_jacobian(test1, v; relstep=1e-3)
end

# ╔═╡ 9b52da32-5fed-4c9a-b98e-c681ad9e9629
isapprox.(g1[:,1], g1[:,3], rtol = 1e-2, atol = 1e-2)  #Forwarddiff

# ╔═╡ 839ff456-c67e-4b52-b95d-673d0df678b9
isapprox.(g1[:,2], g1[:,3], rtol = 1e-2, atol = 1e-2) #Zygote

# ╔═╡ 83fc590c-1586-45b2-b6ed-531459612438
g1[:,1] .- g1[:,3]

# ╔═╡ a326d6e0-5aab-42e0-ae6c-95a83ed9e3c1
g1[:,2] .- g1[:,3]

# ╔═╡ 8c299d82-2b9b-4a4d-9045-eb5ed5ce1b72
g1

# ╔═╡ 446c1f4c-da56-4158-b645-bf38909736f4
begin
	g2 = zeros(length(v),3)
	g2[:,1] = ForwardDiff.gradient(test2, v)
	g2[:,2] = Zygote.gradient(test2, v)[1]
	g2[:,3] = FiniteDiff.finite_difference_jacobian(test2, v; relstep=1e-3)
end

# ╔═╡ 9a7cb749-8634-441d-b5a0-62fb1513a884
isapprox.(g2[:,1], g2[:,3], rtol = 1e-2, atol = 1e-2)  #Forwarddiff

# ╔═╡ 6b02d194-135f-471d-888e-608601fd62cf
isapprox.(g2[:,2], g2[:,3], rtol = 1e-2, atol = 1e-2) #Zygote

# ╔═╡ e7f12bff-ebd3-42e7-b2ec-b27a1b48707e
g2

# ╔═╡ 8156153c-4e97-40c1-8cf7-a3220c68598d
isapprox.(g2[:,1], g2[:,3], rtol = 1e-2, atol = 1e-2)  #Forwarddif

# ╔═╡ fb6fbf19-60dc-46f2-8a4e-493f2988b24c
isapprox.(g2[:,2], g2[:,3], rtol = 1e-2, atol = 1e-2) #Zygote

# ╔═╡ d202a312-0c2f-4c60-964f-87264a8a741d
isapprox.(g2[:,1], g2[:,3], rtol = 1e-2, atol = 1e-2)  #Forwarddif

# ╔═╡ b0984ff2-ecc3-4d34-82bf-54ec9002ee06
isapprox.(g2[:,2], g2[:,3], rtol = 1e-2, atol = 1e-2) #Zygote

# ╔═╡ 6e3169c7-6033-46b7-80c3-4d1a48de7184
isapprox.(g2[:,1], g2[:,3], rtol = 1e-2, atol = 1e-2)  #Forwarddif

# ╔═╡ 13ecf7d8-45f4-48b4-b8cf-0957ef0f93fe
isapprox.(g2[:,2], g2[:,3], rtol = 1e-2, atol = 1e-2) #Zygote

# ╔═╡ efbb1cd1-4077-4dab-abe9-15870fa78f16
begin
	g3 = zeros(length(v),3)
	g3[:,1] = ForwardDiff.gradient(test3, v)
	g3[:,2] = Zygote.gradient(test3, v)[1]
	g3[:,3] = FiniteDiff.finite_difference_jacobian(test3, v; relstep=1e-3)
end

# ╔═╡ 2b5720d8-fc03-45a7-a98d-dc43adc1cb7e
g3

# ╔═╡ d9b63d9e-6d68-4b2c-9328-9a4c9dd10749
begin
	g4 = zeros(length(v),3)
	g4[:,1] = ForwardDiff.gradient(test4, v)
	g4[:,2] = Zygote.gradient(test4, v)[1]
	g4[:,3] = FiniteDiff.finite_difference_jacobian(test4, v; relstep=1e-3)
end

# ╔═╡ 0f51dce4-24bf-4de5-b933-a5a9962b04ec
g4

# ╔═╡ d3957049-52bf-409e-876c-493ffc40a530
begin
	g5 = zeros(length(v),3)
	g5[:,1] = ForwardDiff.gradient(test5, v)
	g5[:,2] = Zygote.gradient(test5, v)[1]
	g5[:,3] = FiniteDiff.finite_difference_jacobian(test5, v; relstep=1e-3)
end

# ╔═╡ 06b035ad-c228-4291-9e27-72325866ef4d
g5

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
GenericLinearAlgebra = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OMEinsum = "ebe7aa44-baf0-506c-a96f-8464559b3922"
StatsFuns = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[compat]
FiniteDiff = "~2.8.1"
ForwardDiff = "~0.10.23"
GenericLinearAlgebra = "~0.2.7"
OMEinsum = "~0.6.2"
StatsFuns = "~0.9.12"
Zygote = "~0.6.30"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "e527b258413e0c6d4f66ade574744c94edef81f8"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.40"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[BFloat16s]]
deps = ["LinearAlgebra", "Printf", "Random", "Test"]
git-tree-sha1 = "a598ecb0d717092b5539dbbe890c98bac842b072"
uuid = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"
version = "0.2.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BatchedRoutines]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8ee75390ba4bbfaf9aa48c121857b0da9a914265"
uuid = "a9ab73d0-e05c-5df1-8fde-d6a4645b8d8e"
version = "0.2.1"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[CUDA]]
deps = ["AbstractFFTs", "Adapt", "BFloat16s", "CEnum", "CompilerSupportLibraries_jll", "ExprTools", "GPUArrays", "GPUCompiler", "LLVM", "LazyArtifacts", "Libdl", "LinearAlgebra", "Logging", "Printf", "Random", "Random123", "RandomNumbers", "Reexport", "Requires", "SparseArrays", "SpecialFunctions", "TimerOutputs"]
git-tree-sha1 = "2c8329f16addffd09e6ca84c556e2185a4933c64"
uuid = "052768ef-5323-5732-b1bb-66c8b64840ba"
version = "3.5.0"

[[ChainRules]]
deps = ["ChainRulesCore", "Compat", "LinearAlgebra", "Random", "RealDot", "Statistics"]
git-tree-sha1 = "035ef8a5382a614b2d8e3091b6fdbb1c2b050e11"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.12.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

[[ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "3287dacf67c3652d3fed09f4c12c187ae4dbb89a"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.4.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "6406b5112809c08b1baa5703ad274e1dded0652f"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.23"

[[GPUArrays]]
deps = ["Adapt", "LinearAlgebra", "Printf", "Random", "Serialization", "Statistics"]
git-tree-sha1 = "7772508f17f1d482fe0df72cabc5b55bec06bbe0"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "8.1.2"

[[GPUCompiler]]
deps = ["ExprTools", "InteractiveUtils", "LLVM", "Libdl", "Logging", "TimerOutputs", "UUIDs"]
git-tree-sha1 = "77d915a0af27d474f0aaf12fcd46c400a552e84c"
uuid = "61eb1bfa-7361-4325-ad38-22787b887f55"
version = "0.13.7"

[[GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "ac44f4f51ffee9ff1ea50bd3fbb5677ea568d33d"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.2.7"

[[IRTools]]
deps = ["InteractiveUtils", "MacroTools", "Test"]
git-tree-sha1 = "95215cd0076a150ef46ff7928892bc341864c73c"
uuid = "7869d1d1-7146-5819-86e3-90919afe41df"
version = "0.4.3"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "46092047ca4edc10720ecab437c42283cd7c44f3"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.6.0"

[[LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6a2af408fe809c4f1a54d2b3f188fdd3698549d6"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.11+0"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OMEinsum]]
deps = ["AbstractTrees", "BatchedRoutines", "CUDA", "ChainRulesCore", "Combinatorics", "LinearAlgebra", "MacroTools", "Requires", "Test", "TupleTools"]
git-tree-sha1 = "84907cc6e34da3f3adb8b6ecfb38cb9f9a75cb56"
uuid = "ebe7aa44-baf0-506c-a96f-8464559b3922"
version = "0.6.2"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Random123]]
deps = ["Libdl", "Random", "RandomNumbers"]
git-tree-sha1 = "0e8b146557ad1c6deb1367655e052276690e71a3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.4.2"

[[RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "e7bc80dc93f50857a5d1e3c8121495852f407e6a"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsFuns]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "95072ef1a22b057b1e80f73c2a89ad238ae4cfff"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.12"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "7cb456f358e8f9d102a8b25e8dfedf58fa5689bc"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.13"

[[TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zygote]]
deps = ["AbstractFFTs", "ChainRules", "ChainRulesCore", "DiffRules", "Distributed", "FillArrays", "ForwardDiff", "IRTools", "InteractiveUtils", "LinearAlgebra", "MacroTools", "NaNMath", "Random", "Requires", "SpecialFunctions", "Statistics", "ZygoteRules"]
git-tree-sha1 = "2c30f2df0ba43c17e88c8b55b5b22c401f7cde4e"
uuid = "e88e6eb3-aa80-5325-afca-941959d7151f"
version = "0.6.30"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─fc634506-632f-47ba-8a07-b5ce108c2012
# ╟─40e6e0a9-ca69-47ba-b982-ea2e9d098578
# ╟─56a94c71-1329-40cd-8bed-3748307ba2b4
# ╟─fe3764e3-f309-41f5-8d2d-8d1783d87cba
# ╟─604dfba0-5d9d-4fbf-82bc-7f19f4032c5b
# ╟─9b52da32-5fed-4c9a-b98e-c681ad9e9629
# ╠═839ff456-c67e-4b52-b95d-673d0df678b9
# ╠═83fc590c-1586-45b2-b6ed-531459612438
# ╠═a326d6e0-5aab-42e0-ae6c-95a83ed9e3c1
# ╠═8c299d82-2b9b-4a4d-9045-eb5ed5ce1b72
# ╟─22066571-3e2b-43a1-862c-6a40bc3f4b47
# ╟─d5b61aa8-5125-4a3d-804f-26b55879c0f5
# ╠═97e58e64-7dfc-4414-a7e7-8d1f38e9a4ec
# ╠═446c1f4c-da56-4158-b645-bf38909736f4
# ╠═9a7cb749-8634-441d-b5a0-62fb1513a884
# ╠═6b02d194-135f-471d-888e-608601fd62cf
# ╠═e7f12bff-ebd3-42e7-b2ec-b27a1b48707e
# ╟─336ee5a5-eee2-416d-a72e-27cdff89369c
# ╠═ad7aa9c3-d546-4e1e-931a-d56377532cf2
# ╠═10516ffa-6e74-4d95-9ebd-492166e8ca28
# ╠═efbb1cd1-4077-4dab-abe9-15870fa78f16
# ╠═8156153c-4e97-40c1-8cf7-a3220c68598d
# ╠═fb6fbf19-60dc-46f2-8a4e-493f2988b24c
# ╠═2b5720d8-fc03-45a7-a98d-dc43adc1cb7e
# ╟─36dea9c1-f358-4bfc-b243-2c33919ae5c6
# ╠═123622a5-f727-4426-9c6e-5fdf749ec610
# ╠═0208996b-a133-4c72-9a5f-048f963924ac
# ╠═d9b63d9e-6d68-4b2c-9328-9a4c9dd10749
# ╠═d202a312-0c2f-4c60-964f-87264a8a741d
# ╠═b0984ff2-ecc3-4d34-82bf-54ec9002ee06
# ╠═0f51dce4-24bf-4de5-b933-a5a9962b04ec
# ╟─9bc2dd49-ee54-4e58-bac0-13426d9a7098
# ╠═f50f2e48-3da1-49e5-9f2e-c3eb22e3f6b4
# ╠═23065c11-ebe1-47a2-a770-c4283dd2c36b
# ╟─d3957049-52bf-409e-876c-493ffc40a530
# ╠═6e3169c7-6033-46b7-80c3-4d1a48de7184
# ╠═13ecf7d8-45f4-48b4-b8cf-0957ef0f93fe
# ╠═06b035ad-c228-4291-9e27-72325866ef4d
# ╠═4de5fbf3-afce-4d12-883d-ab3b3b6da1b6
# ╟─88ae1058-abcf-4567-86dd-8e1f1c5f7871
# ╟─2bbca5ed-12bc-4646-a401-d691e9c138eb
# ╟─91fc3a50-0420-4e3f-8e15-58fe446ce0be
# ╠═769b2024-6459-44c2-ab36-01b74944b921
# ╠═f3e91cc8-41cc-11ec-0cc4-0fe7d3ecfd69
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
