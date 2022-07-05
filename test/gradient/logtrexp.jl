#logtrexp
for solver in [cpu_solver, gpu_solver] 
        D = rand(2:6)
        M = rand(D, D) |> symmetrize
        t = rand(1)
        ∂M = x->logtrexp(t[1], x)

        gradt_ndiff = ngradient(x->logtrexp(x...,M), t)[1] 
        gradM_ndiff = ngradient(x->logtrexp(t[1], x), M)[1] 

    @testset "$(solver)" begin 
        gradt_zygote = solver(m->Zygote.gradient(x->logtrexp(x...,m), t), M)[1] 
        gradM_zygote = solver(m->Zygote.gradient(x->logtrexp(t[1], x), m), M)[1] |> Array
        @test ≈(gradt_zygote, gradt_ndiff; rtol = 1e-5, atol = 1e-5)
        @test ≈(gradM_zygote, gradM_ndiff; rtol = 1e-5, atol = 1e-5)
    end
end