#free energy
for solver in [cpu_solver, gpu_solver] 
    @testset "$(solver)" begin
        @testset "TFIsing model" begin
            β = 2.0; χ = 2
            m = TFIsing(1.0, 1.0)
            ψ = init_cmps(χ) |> diagQ
            dQ = convert(Vector, diag(ψ.Q))
            R = convert(Array, ψ.R)
            
            Tm = solver(x->x, m.Tmatrix)
            loss() = free_energy(solver(CMPS_generate, diagm(dQ), R), Tm, β)
            pars = Zygote.Params([dQ, R])
            p0, f, g! = optim_functions(loss, pars)
                
            grad_ndiff = ngradient(f, p0)[1]
            grad_zygote = similar(p0); g!(grad_zygote, p0)
            @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
        end
    end
end