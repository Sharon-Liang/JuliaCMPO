using JuliaCMPO, Test
using CUDA; CUDA.allowscalar(false)
processor = GPU
solver = solver_function(processor)

@show solver

@testset "otimes.jl" begin
    include("otimes.jl")
end

@testset "CTensorProducts.jl" begin
    include("CTensorProducts.jl")
end

@testset "CMPSMatrix.jl" begin
    include("CMPSMatrix.jl")
end

@testset "logtrexp.jl" begin
    include("logtrexp.jl")
end

@testset "CMPSOperations.jl" begin
    include("CMPSOperations.jl")
end

@testset "ModelConstruct.jl" begin
    @testset "Ut' * T * Ut = transpose(T)" begin
        wid = 3
        for m in [TFIsing(1.0,1.0), XYmodel(), XYmodel_2D_helical(1, expand=true), XXZmodel_2D_helical(2.0, wid, expand=true)]
            Tm = solver(x->x, m.Tmatrix)
            Ut = solver(x->x, m.Ut)
            @test  Ut' * Tm * Ut == transpose(Tm)
         end
    end
end

#gradient test
@testset "gradient/logtrexp.jl" begin
    include("./gradient/logtrexp.jl")
end

@testset "gradient/logfidelity.jl" begin
    include("./gradient/logfidelity.jl")
end

@testset "gradient/project.jl" begin
    include("./gradient/project.jl")
end

@testset "gradient/free_energy.jl" begin
    include("gradient/free_energy.jl")
end

@testset "CMPSInitiate.jl" begin
    include("CMPSInitiate.jl")
end