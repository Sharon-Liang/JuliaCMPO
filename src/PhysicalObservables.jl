#module PhysicalObservables
#include("Setup.jl")

"""
NN Transvers field Ising model
    H = ∑ J Zi Zj + ∑ Γ Xi
"""
function TFIsing(J::Real, Γ::Real; field='N')
    if field == 'N'
        h = zeros(2,2)
    else
        η = 1.e-5
        h = η .* pauli(field)
    end
    return cmpo(Γ*pauli('x')+h, √J*pauli('z'), √J*pauli('z'), zeros(2,2))
end

"""
Free energy
"""
function free_energy(ψ::cmps, W::cmpo, β::Real)
    K = ψ * W * ψ ; H = ψ * ψ
    res = logtrexp(-β*K)- logtrexp(-β*H)
    return -1/β * res
end

function free_energy(param::Array{Float64,3}, W::cmpo, β::Real)
    ψ = cmps(param[:,:,1], param[:,:,2])
    free_energy(ψ, W, β)
end


"""
The thermal average of local opeartors
"""
function thermal_average(Op::AbstractArray, ψ::cmps, W::cmpo, β::Real)
    eye = Matrix(1.0I, size(ψ.Q))
    Op = eye ⊗ Op ⊗ eye
    K = ψ * W * ψ |> symmetrize |> Hermitian
    e, v = eigen(-β*K)
    m = maximum(e)
    Op = v' * Op * v
    den = exp.(e .- m) |> sum
    num = exp.(e .- m) .* diag(Op) |> sum
    return num/den
end

"""
The local two-time correlation functions
"""
function correlation_2time(A::AbstractArray,B::AbstractArray, τ::Number,
                            ψ::cmps, W::cmpo, β::Real)
       eye = Matrix(1.0I, size(ψ.Q))
       A = eye ⊗ A ⊗ eye
       B = eye ⊗ B ⊗ eye
       K = ψ * W * ψ |> symmetrize |> Hermitian
       e, v = eigen(-β * K)
       m = maximum(val)
       A = v' * A * v
       B = v' * B * v
       den = exp.(e .- m) |> sum
       num = 0.0
       for i = 1: length(e), j = 1: length(e)
           num += exp(e[i]-m + τ*(e[i] - e[j])) * A[i,j] * B[j,i]
       end
       return num/den
end

function spectrum(A::AbstractArray,B::AbstractArray, ω::Real,
                    ψ::cmps, W::cmpo, β::Real; type = "ret")
       η = 1.e-5
       if type=="ret" ω = ω + η * 1im
       elseif type=="adv" ω = ω - η * 1im
       elseif type=="ord" ω = ω - η * 1im * sign(ω)
       elseif type=="img" ω = ω * 1im
       else error("type should be 'ret','adv','ord','img'")
       end
       eye = Matrix(1.0I, size(ψ.Q))
       A = eye ⊗ A ⊗ eye
       B = eye ⊗ B ⊗ eye
       K = ψ * W * ψ |> symmetrize |> Hermitian
       e, v = eigen(-β * K)
       m = maximum(val)
       A = v' * A * v
       B = v' * B * v
       den = exp.(e .- m) |> sum
       num = 0.0
       for i = 1: length(e), j = 1: length(e)
           up = exp(e[j]-m) - exp(e[i]-m)
           up = up * A[i,j] * B[j,i]
           down = ω + e[i] - e[j]
           num += up/down
       end
       return num/den
end
#end  # module PhysicalObservables
