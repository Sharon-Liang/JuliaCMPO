#project function
for solver in [cpu_solver, gpu_solver] 
    @testset "$(solver)" begin
        β = rand()
        χ = 4
        D = 2
        loss(v, ψ) = logfidelity(project(ψ, v), ψ, β)
            
        for vir_dim = 1:2
            u = solver(x->x, rand(χ,D))
            ψ0 = solver(x->x, init_cmps(χ, vir_dim))
            grad_ndiff = ngradient(u->loss(Array(u), CTensor(ψ0)), Array(u))[1]  
            grad_zygote = Zygote.gradient(u->loss(u, ψ0), u)[1] |> Array
            @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
        end
    end
end