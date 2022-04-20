using cMPO, Test
using Random; Random.seed!()
using Zygote
using LinearAlgebra

function ngradient(f, xs::AbstractArray...)
    #https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl
    grads = zero.(xs)
    for (x, Δ) in zip(xs, grads), i in 1:length(x)
      δ = sqrt(eps())
      tmp = x[i]
      x[i] = tmp - δ/2
      y1 = f(xs...)
      x[i] = tmp + δ/2
      y2 = f(xs...)
      x[i] = tmp
      Δ[i] = (y2-y1)/δ
    end
    return grads
end

for solver in [cpu_solver, gpu_solver] 
    @testset "$(solver)" begin
        @testset "logtrexp" begin
            D = rand(2:6)
            A = rand(D, D)
            grad_ndiff = ngradient(logtrexp, A)[1] 
            grad_zygote = solver(x->Zygote.gradient(logtrexp, x), A)[1] |> Array
            @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
        end

        @testset "logfidelity" begin
            β = 2.0
            χ = 4
            for vir_dim = 1:2
                ψ0 = solver(x->x, init_cmps(χ+2, vir_dim) )

                ψ = init_cmps(χ, vir_dim) |> diagQ
                dQ = convert(Vector, diag(ψ.Q))
                R = convert(Array, ψ.R)

                loss() = -logfidelity(solver(CMPS_generate, consist_diagm(dQ), R), ψ0, β)
                pars = Zygote.Params([dQ, R])
                p0, f, g! = optim_functions(loss, pars)

                grad_ndiff = ngradient(f, p0)[1]
                grad_zygote = similar(p0); g!(grad_zygote, p0)
                @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
            end
        end


        @testset "project function" begin
            β = 2.0
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

        @testset "free energy: TFIsing model" begin
            β = 2.0; χ = 2
            m = TFIsing(1.0, 1.0)
            ψ = init_cmps(χ) |> diagQ
            dQ = convert(Vector, diag(ψ.Q))
            R = convert(Array, ψ.R)
            
            Tm = solver(x->x, m.Tmatrix)
            loss() = free_energy(solver(CMPS_generate, consist_diagm(dQ), R), Tm, β)
            pars = Zygote.Params([dQ, R])
            p0, f, g! = optim_functions(loss, pars)
                
            grad_ndiff = ngradient(f, p0)[1]
            grad_zygote = similar(p0); g!(grad_zygote, p0)
            @test ≈(grad_zygote, grad_ndiff; rtol = 1e-5, atol = 1e-5)
        end

    end
end












