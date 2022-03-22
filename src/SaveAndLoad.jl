#module SaveAndLoad
function saveCMPS(path::String, ψ::CMPS, dict::Dict = Dict())
    h5open(path, "w") do file
        file["Q"] = ψ.Q
        file["R"] = ψ.R
        for key in keys(dict) 
            file[key] = dict[key] 
        end
    end
end

function readCMPS(path::String; python::Bool = false)
    h5open(path, "r") do file
        Q = read(file["Q"]) 
        R = read(file["R"])
        if python
            Q = Q |> diagm
            R = ein"ijk->jik"(R)
        end
        return CMPS(Q, R)
    end
end


#Reference: https://github.com/JuliaNLSolvers/Optim.jl/blob/master/src/types.jl
"""
    Write Optim.optimize result
"""
function Base.write(io::IO, r::Optim.MultivariateOptimizationResults)
    take = Iterators.take
    failure_string = "failure"
    if Optim.iteration_limit_reached(r)
        failure_string *= " (reached maximum number of iterations)"
    end
    if Optim.f_increased(r) && !Optim.iteration_limit_reached(r)
        failure_string *= " (objective increased between iterations)"
    end
    if isa(r.ls_success, Bool) && !r.ls_success
        failure_string *= " (line search failed)"
    end
    if Optim.time_run(r) > Optim.time_limit(r)
        failure_string *= " (exceeded time limit of $(Optim.time_limit(r)))"
    end

    @printf io " * Status: %s\n\n" Optim.converged(r) ? "success" : failure_string

    @printf io " * Candidate solution\n"
    @printf io "    Final objective value:     %e\n" Optim.minimum(r)
    @printf io "\n"

    @printf io " * Found with\n"
    @printf io "    Algorithm:     %s\n" Optim.summary(r)


    @printf io "\n"
    @printf io " * Convergence measures\n"
    if isa(r.method, NelderMead)
        @printf io "    √(Σ(yᵢ-ȳ)²)/n %s %.1e\n" Optim.g_converged(r) ? "≤" : "≰" Optim.g_tol(r)
    else
        @printf io "    |x - x'|               = %.2e %s %.1e\n"  Optim.x_abschange(r) Optim.x_abschange(r)<=Optim.x_abstol(r) ? "≤" : "≰" Optim.x_abstol(r)
        @printf io "    |x - x'|/|x'|          = %.2e %s %.1e\n"  Optim.x_relchange(r) Optim.x_relchange(r)<=Optim.x_reltol(r) ? "≤" : "≰" Optim.x_reltol(r)
        @printf io "    |f(x) - f(x')|         = %.2e %s %.1e\n"  Optim.f_abschange(r) Optim.f_abschange(r)<=Optim.f_abstol(r) ? "≤" : "≰" Optim.f_abstol(r)
        @printf io "    |f(x) - f(x')|/|f(x')| = %.2e %s %.1e\n"  Optim.f_relchange(r) Optim.f_relchange(r)<=Optim.f_reltol(r) ? "≤" : "≰" Optim.f_reltol(r)
        @printf io "    |g(x)|                 = %.2e %s %.1e\n" Optim.g_residual(r) Optim.g_residual(r)<=Optim.g_tol(r) ?  "≤" : "≰" Optim.g_tol(r)
    end

    @printf io "\n"

    @printf io " * Work counters\n"
    @printf io "    Seconds run:   %d  (vs limit %d)\n" Optim.time_run(r) isnan(Optim.time_limit(r)) ? Inf : Optim.time_limit(r)
    @printf io "    Iterations:    %d\n" Optim.iterations(r)
    @printf io "    f(x) calls:    %d\n" Optim.f_calls(r)
    if !(isa(r.method, NelderMead) || isa(r.method, SimulatedAnnealing))
        @printf io "    ∇f(x) calls:   %d\n" Optim.g_calls(r)
    end
    if isa(r.method, Newton) || isa(r.method, NewtonTrustRegion)
        @printf io "    ∇²f(x) calls:  %d\n" Optim.h_calls(r)
    end
    return
end


"""
    write Optim trace
"""
function Base.write(io::IO, tr::Optim.OptimizationTrace)
    @printf io "\n"
    @printf io "Iter      Function value     Gradient norm  \n"
    @printf io "------   ----------------   ----------------\n"
    for state in tr
        writedlm(io, [state.iteration state.value state.g_norm])
    end
    return
end



"""
    write MERA update trace
"""
function Base.write(io::IO, tr::MeraUpdateTrace)
    @printf io "step          θ            loss_diff           fidelity    \n"
    @printf io "----  ----------------  ----------------   ----------------\n"
    for state in tr
        writedlm(io, [state.SN state.θ state.loss_diff state.fidelity])
    end
    return
end

function Base.println(step::MeraUpdateStep)
    str = @sprintf "%3i   %.16f   %.16e   %.16f \n" step.SN step.θ step.loss_diff step.fidelity
    println(str)
    return
end

#end    # module SaveAndLoad