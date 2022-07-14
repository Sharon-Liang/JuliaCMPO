#logfidelity
for solver in [cpu_solver, gpu_solver] 
    @testset "$(solver)" begin
        β = 2.0
        χ = 4
        for vir_dim = 1:2
            ψ0 = solver(x->x, init_cmps(χ+2, vir_dim) )

            ψ = init_cmps(χ, vir_dim) |> diagQ
            dQ = convert(Vector, diag(ψ.Q))
            R = convert(Array, ψ.R)

            loss() = -logfidelity(solver(CMPS_generate, diagm(dQ), R), ψ0, β)
            pars = Zygote.Params([dQ, R])
            p0, f, g! = optim_functions(loss, pars)

            grad_ndiff = ngradient(f, p0)[1]
            grad_zygote = similar(p0); g!(grad_zygote, p0)
            @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
        end
    end
end