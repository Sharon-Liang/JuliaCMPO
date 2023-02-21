#=
### *Core* : *MERA update*
=#
"""
    MeraUpdateOptions

Data structure of keyword arguments of mera_update function.

Keyword Arguments:
* atol : Tolerance of absolute difference of the fidelity between two normalized cMPS and `1.0`, defalt `atol=1.e-5`.
* btol : Tolerance of absolute difference of logfidelity between two MERA update steps, defalt `btol = 1.e-12`.
* maxiter : Maximum iteration step, defalt 50.
* interpolate : `Bool`, whether to interpolate between two MERA update steps, defalt `true`.

See also: [`mera_update`](@ref)
"""
@with_kw struct MeraUpdateOptions
    atol::Float64 = 1.e-5
    btol::Float64 = 1.e-12
    maxiter::Int64 = 50
    interpolate::Bool = true
end


"""
    _interpolate_isometry(p1, p2, θ)

Interpolate between two unitrary matrices. One gets `p1` if ``θ = π/2`` and gets `p2` if ``θ = 0``.

Reference: 
* https://mathoverflow.net/questions/262560/natural-ways-of-interpolating-unitary-matrices
* https://groups.google.com/forum/#!topic/manopttoolbox/2zhx67doXaU
"""
function _interpolate_isometry(p1::AbstractMatrix, p2::AbstractMatrix, θ::Real)
    mix = sin(θ) * p1 + cos(θ) * p2
    F = svd(mix)
    return F.U * F.Vt 
end



"""
    mera_update(ψ₀::CMPS, χ::Integer, β::Real[; atol::Float64=1.e-5, btol::Float64=1.e-12, maxiter::Int=50, interpolate::Bool=true])

To Compress a cMPS `ψ₀` to bond dimension `χ` by inserting isometries. And using iterative SVD update with line search to find the optimal isometry.

Arguments:
* ψ₀ : The cMPS to be compressed.
* χ  : Target bond dimension.
* β  : Inverse temperature.

Keyword Arguments:
* atol : Tolerance of absolute difference of the fidelity between two normalized cMPS and `1.0`, defalt `atol=1.e-5`.
* btol : Tolerance of absolute difference of logfidelity between two MERA update steps, defalt `btol = 1.e-12`.
* maxiter : Maximum iteration step, defalt 50.
* interpolate : `Bool`, whether to interpolate between two MERA update steps, defalt `true`.

See also: [`MeraUpdateOptions`](@ref)
"""
function mera_update(ψ₀::CMPS, χ::Integer, β::Real; atol::Float64=1.e-5, btol::Float64=1.e-12, maxiter::Int=50, interpolate::Bool=true)
    #count iteration step
    step = 1
    n₀ = norm(ψ₀, β)

    #Define loss function: P is the isometry
    loss(P) = logfidelity(project(ψ₀, P), ψ₀, β, false)

    #Calculate the current isometry Pc corresponding to current ψ₀
    Q = symmetrize(ψ₀.Q)
    _, v = eigensolver(Q)
    Pc = v[:, end-χ+1:end]
    
    #Initiate loss function value
    Lp = 9.9e9
    Lc = loss(Pc)

    #Calculate the difference between current fidelity and 1.0
    ΔF = abs(exp(Lc)/n₀ - 1.0)
    #Calculate the difference of logfidelity between the current step and the previous step
    ΔlnF = abs(Lc - Lp)
    
    println("Adaptive MERA Update\n")
    println("step        θ/π               ΔlnF                1.0 - F      ")     
    println("----  ----------------  -----------------   -------------------")
    println(@sprintf "%03i   %.10e   %.10e   %.10e" step 1.0 ΔlnF ΔF)

    while step < maxiter
        step += 1   
        grad = Zygote.gradient(loss, Pc)[1]
        F = svd(grad)
        Pn = F.U * F.Vt
 
        #interpolate between unitary matrices
        θ = π
        proceed = interpolate
        while proceed
            θ = θ/2
            #12-times bisection, cos(θ) = 0.9999989926433588
            if θ < π/(1.9^12) 
                Pn = Pc
                proceed = false
            else
                Pi = _interpolate_isometry(Pn, Pc, θ)
                Li = loss(Pi)
                if Li > Lc
                    Pn = Pi
                    Lc = Li
                    proceed = false
                end
            end     
        end

        ΔF = abs(exp(Lc)/n₀ - 1.0)
        ΔlnF = abs(Lc - Lp)
        println(@sprintf "%03i   %.10e   %.10e   %.10e" step θ/π ΔlnF ΔF)

        #Update
        Pc = Pn
        Lp = Lc

        if ΔF < atol || ΔlnF < btol break end
    end
    return project(ψ₀, Pc)
end

function mera_update(ψ₀::CMPS, χ::Integer, β::Real, options::MeraUpdateOptions)
    @unpack atol, btol, maxiter, interpolate = options
    return mera_update(ψ₀, χ, β; atol, btol, maxiter, interpolate)
end




#=
### *Core* : *Compress*
=#
"""
    CompressOptions

Data structure of keyword arguments of compress_cmps function.

Keyword Arguments:
* processor : Precessor used, defalt `CPU`.
* mera_update_options : Options used in mera_update function.
* optim_options : Optins used in optimization process.

See also: [`compress_cmps`](@ref)
"""
@with_kw struct CompressOptions
    processor::Processor = CPU
    mera_update_options::MeraUpdateOptions = MeraUpdateOptions()
    optim_options::Optim.Options = Optim.Options(
        f_tol = 2.220446049250313e-9, 
        g_tol = 1.e-5, 
        iterations = 100, 
        show_trace = true)
end


"""
    compress_cmps(ψ₀, χ, β[, init::Union{Nothing, CMPS} = nothing]; options::CompressOptions)

Compress a CMPS `ψ₀` to bond dimension `χ` at temperature `β`.

Arguments:
* ψ₀ : The cMPS to be compressed.
* χ  : Target bond dimension.
* β  : Inverse temperature.
"""
function compress_cmps(ψ₀::CMPS, χ::Integer, β::Real, init::Union{Nothing, CMPS} = nothing; compress_options::CompressOptions = CompressOptions())
    @unpack processor, mera_update_options, optim_options = compress_options

    #Generate solver function and initiate cMPS accordingly
    ψ₀ = solver(ψ₀)

    init === nothing ? ψ = mera_update(ψ₀, χ, β, mera_update_options) : ψ = solver(init)
    @assert size(ψ.Q) == (χ, χ) 

    #Calculate initial fidelity
    Fi = fidelity(ψ, ψ₀, β, true)

    ψ = ψ |> diagQ
    #This conversion is necessary because parameters in type Array is required in Optim.optimize function.  
    Qd = convert(Vector, diag(ψ.Q))
    R  = convert(Array, ψ.R)

    #Generate loss function and its gradient function
    loss() = -logfidelity(solver(CMPS(diagm(Qd), R)), ψ₀, β, false)
    p0, f, g! = optim_functions(loss, Params([Qd, R]))
    
    #Optimize
    res = Optim.optimize(f, g!, p0, LBFGS(), optim_options)
    
    ψ = solver(CMPS(diagm(Qd), R))
    #calculate final fidelity
    Ff = fidelity(ψ, ψ₀, β, true)

    println(res)
    println(@sprintf "|1 - Fidelity| Change: %.5e -> %.5e\n" 1-Fi 1-Ff)

    return ψ
end


#=
### *Core* : *Initiate*
=#
"""
    init_cmps(χ[, D = 1, hermitian::Bool=true])

Randomly initiate a cMPS local tensor with bond dimension `D+1` and virtual bond dimension `χ`.
"""
function init_cmps(χ::Int64, D::Int64 = 1, hermitian::Bool = true)
    Q = rand(χ, χ)
    D == 1 ? R = rand(χ, χ) : R = rand(χ, χ, D)

    if hermitian
        Q = symmetrize(Q)
        for d = 1:D
            R[:,:,d] = symmetrize(R[:,:,d])
        end
    end
    return CMPS(Q, R)
end


"""
    init_cmps(χ::Integer, Tₘ::CMPO, β::Real; options::CompressOptions)

Using the boundary of a cMPO to initiate a cMPS with bond dimension `χ` at inverse temperature `β`.
"""
function init_cmps(χ::Integer, Tₘ::CMPO, β::Real; compress_options::CompressOptions = CompressOptions())
    ψ = CMPS(Tₘ.Q, Tₘ.R)
    while size(ψ.Q, 1) < χ
        ψ = Tₘ * ψ 
    end

    if size(ψ.Q, 1) > χ 
        ψ = compress_cmps(ψ, χ, β; compress_options)
    end

    return ψ
end