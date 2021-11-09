import ChainRules

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