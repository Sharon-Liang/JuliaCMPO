#import ChainRules
#import OMEinsum:einsum_grad, _insertat, DynamicEinCode
"""
function ChainRules.rrule(::typeof(GenericLinearAlgebra.eigvals), A; kwargs...)
    F, eigen_back = ChainRules.rrule(GenericLinearAlgebra.eigen, A; kwargs...)
    λ = F.values
    function eigvals_pullback(Δλ)
        ∂F = ChainRules.Tangent{typeof(F)}(values = Δλ)
        _, ∂A = eigen_back(∂F)
        return ChainRules.NoTangent(), ∂A
    end
    return λ, eigvals_pullback
end

function ChainRules.rrule(
    ::typeof(GenericLinearAlgebra.eigen),
    A;
    kwargs...,
    )
    F = GenericLinearAlgebra.eigen(A; kwargs...)
    function eigen_pullback(ΔF::ChainRules.Tangent)
        λ, U = F.values, F.vectors
        Δλ, ΔU = ΔF.values, ΔF.vectors
        ΔU = ΔU isa ChainRules.AbstractZero ? ΔU : copy(ΔU)
        ∂A = ChainRules.eigen_rev!(A, λ, U, Δλ, ΔU)
        return ChainRules.NoTangent(), ∂A
    end
    eigen_pullback(ΔF::ChainRules.AbstractZero) = (ChainRules.NoTangent(), ΔF)
    return F, eigen_pullback
end

function einsum_grad(ixs, @nospecialize(xs), iy, size_dict, cdy, i)
    nixs = _insertat(ixs, i, iy)
    nxs  = _insertat( xs, i, cdy)
    niy = ixs[i]
    y = einsum(DynamicEinCode(nixs, niy), nxs, size_dict)
    y = conj(y)  # do not use `conj!` to help computing Hessians.
    typeof(y) == typeof(xs[i]) && return y
    #xs[i] isa Array{<:Real} && return convert(typeof(xs[i]), real(y))
    #convert(typeof(xs[i]), y)
end
"""

"""
    rrule for logtrexp function, eltype(M) <: Real
"""
function ChainRules.rrule(::typeof(logtrexp), 
            M::AbstractArray; 
            device::Symbol=:cpu)
    e, v = symeigen(M, device = device)
    y = logsumexp(e)
    function logtrexp_pullback(ȳ)
        ∂y_∂M = v * diagm(exp.(e .- y)) * v' |> transpose
        M̄ = ȳ * ∂y_∂M
        return ChainRules.NoTangent(), M̄, ChainRules.ZeroTangent()
    end
    return y, logtrexp_pullback
end



