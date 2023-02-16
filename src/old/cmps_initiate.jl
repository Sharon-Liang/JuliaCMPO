#module cMPSCompress
"""
    init_cmps([T::DataType = Float64,] χ[, vD = 1; hermitian::Bool=true])

Randomly initiate a `CMPS` with bond dimension `χ` and virtual bond dimension `vD`.
"""
function init_cmps(dtype::DataType, χ::Int64, vD::Int64 = 1; hermitian::Bool = true)
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
    return cmps_generate(Q,R)
end
init_cmps(χ::Int64, vD::Int64; kwargs...) = init_cmps(Float64, χ, vD; kwargs...)
init_cmps(χ::Int64; kwargs...) = init_cmps(Float64, χ, 1; kwargs...)


"""
    init_cmps(χ::Integer, Tm::AbstractCMPO, β::Real; options::CompressOptions)

Initiate a `CMPS` accoring to the `CMPO` tensor `Tm`.
"""
function init_cmps(χ::Integer, Tm::AbstractCMPO, β::Real; compress_options::CompressOptions)
    if compress_options.show_trace println(_header("Initiate CMPS")) end
    ψ = cmps_generate(Tm.Q, Tm.R)
    while size(ψ.Q, 1) < χ
        ψ = Tm * ψ 
    end

    if size(ψ.Q,1) > χ 
        res = compress_cmps(ψ, χ, β; compress_options)
        ψ = res.cmps
    end
    return ψ
end

#end    # module cMPSCompress