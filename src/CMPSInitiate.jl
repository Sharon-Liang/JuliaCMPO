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
    return CMPS_generate(Q,R)
end


"""
    Initiate via boundary cMPS
"""
function init_cmps(bondD::Integer, Tm::AbstractCMPO, β::Real; 
                   options::CompressOptions=CompressOptions())
    if options.show_trace
        println("----------------------------Initiate CMPS-----------------------------") 
    end
    ψ = CMPS_generate(Tm.Q, Tm.R)
    while size(ψ.Q, 1) < bondD
        ψ = Tm * ψ 
    end

    if size(ψ.Q, 1) > bondD 
        res = compress_cmps(ψ, bondD, β, options=options)
        ψ = res.ψ
    end
    return ψ
end

#end    # module cMPSCompress