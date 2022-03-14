#module cMPSCompress
"""
    Randomly initiate a cMPS with bond dimension χ and virtual bond dimension vD
"""
function init_cmps(χ::Int64, vD::Int64 = 1; 
    hermitian::Bool = true, dtype::DataType = Float64)
    Q = rand(dtype, χ, χ)
    if vD == 1
        R = rand(dtype, χ, χ)
    else
        R = rand(dtype, χ, χ, vD)
    end

    if hermitian
        Q = Q |> symmetrize
        for d = 1:vD
            R[:,:,d] = R[:,:,d] |> symmetrize
        end
    end
    return CMPS(Q,R)
end


"""
    Fidelity between the target cMPS |ψ⟩ and the origional cMPS T|r⟩:
        Fd = ⟨ψ|T|r⟩/√(⟨ψ|ψ⟩)
        fidelity(ψ, ψ0) = ln(Fd)
"""
fidelity(ψ::CMPS, ψ0::CMPS, β::Real) = log_overlap(ψ, ψ0, β) - 0.5*log_overlap(ψ, ψ, β)


"""
    compress a CMPS(ψ0) to a target dimension(χ)
"""
function cmps_compress(ψ0::CMPS, χ::Integer, β::Real; init::Union{CMPS, Nothing} = nothing)
    χ0 = size(ψ0.Q)[1]
    length(size(ψ0.R)) == 2 ? vir_dim = 1 : vir_dim = size(ψ0.R)[3]
    if χ0 <= χ
        @warn "The bond dimension of the initial CMPS ≤ target bond dimension."
        return ψ
    else
        init === nothing ? ψ = init_cmps(χ, vir_dim) : ψ = init
        if size(ψ.Q) != (χ, χ) 
            msg = "χ of init cmps should be $(χ) instead of $(size(ψ.Q)[1])"
            @error DimensionMismatch(msg)
        else
            ψ = ψ |> diagQ
            loss() = -fidelity(CMPS(diag(ψ.Q)|> diagm, ψ.R), ψ0, β)

            #FluxOptTools: Zygote.refresh() currently needed when defining new adjoints
            Zygote.refresh() 
            lossfun, gradfun, fg!, p0 = optfuns(loss, Params([ψ.Q, ψ.R]))
            opt = Optim.optimize(Optim.only_fg!(fg!), p0, LBFGS(),Optim.Options(iterations=10000, store_trace=true))
            return ψ
        end
    end 
end

#end    # module cMPSCompress