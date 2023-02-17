using Test, JuliaCMPO
using LinearAlgebra

@show solver = solver_function(processor)
#=
### *logtrexp* 
=#
@testset "matrix logtrexp" begin
    D = 4; T = ComplexF64
    t = rand()
    for A in [randn(D,D), randn(T,D,D)]
        A = symmetrize(A)
        log_trexp = exp(t*A) |> tr |> real |> log
        @test logtrexp(t, solver(A)) ≈ log_trexp
    end
end