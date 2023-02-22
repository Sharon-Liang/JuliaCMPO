using Test, JuliaCMPO
using JLD, DelimitedFiles
import JuliaCMPO._j1_j2_block

processor = CPU
obsv_functions = [free_energy, energy, entropy]

#=
### *J₁-J₂ Chain* : *model cMPO*
=#
@testset " J₁-J₂ Chain cMPO" begin
    J1 = 1.0; J2 = 0.5
    @testset "PX block" begin
        Tx = _j1_j2_block(J1, J2, PX)
        @test Tx.Q == [0.0 0.0; 0.0 0.0]
        @test Tx.R == [0.0 0.5; 0.5 0.0;;; 0.0 0.0; 0.0 0.0]
        @test Tx.L == [-0.0 -0.5; -0.5 -0.0;;; -0.0 -0.3535533905932738; -0.3535533905932738 -0.0]
        @test Tx.P == [0.0 0.0; 0.0 0.0;;; 0.7071067811865476 0.0; 0.0 0.7071067811865476;;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
    end

    @testset "iPY block" begin
        Ty = _j1_j2_block(-J1, -J2, iPY)
        @test Ty.Q == [0.0 0.0; 0.0 0.0]
        @test Ty.R == [0.0 0.5; -0.5 0.0;;; 0.0 0.0; 0.0 0.0]
        @test Ty.L == [0.0 0.5; -0.5 0.0;;; 0.0 0.3535533905932738; -0.3535533905932738 0.0]
        @test Ty.P == [0.0 0.0; 0.0 0.0;;; 0.7071067811865476 0.0; 0.0 0.7071067811865476;;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
    end

    @testset "PX block" begin
        Tz = _j1_j2_block(J1, J2, PZ)
        @test Tz.Q == [0.0 0.0; 0.0 0.0]
        @test Tz.R == [0.5 0.0; 0.0 -0.5;;; 0.0 0.0; 0.0 0.0]
        @test Tz.L == [-0.5 -0.0; -0.0 0.5;;; -0.3535533905932738 -0.0; -0.0 0.3535533905932738]
        @test Tz.P == [0.0 0.0; 0.0 0.0;;; 0.7071067811865476 0.0; 0.0 0.7071067811865476;;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
    end

    @testset "model cMPO" begin
        Tₘ = model(J1J2Chain(), J1, J2)
        @test Tₘ.Q == [0.0 0.0; 0.0 -0.0]
        @test Tₘ.R == [0.0 0.5; 0.5 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.5; -0.5 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.5 0.0; 0.0 -0.5;;; 0.0 0.0; 0.0 0.0]
        @test Tₘ.L == [-0.0 -0.5; -0.5 -0.0;;; -0.0 -0.3535533905932738; -0.3535533905932738 -0.0;;; 0.0 0.5; -0.5 0.0;;; 0.0 0.3535533905932738; -0.3535533905932738 0.0;;; -0.5 -0.0; -0.0 0.5;;; -0.3535533905932738 -0.0; -0.0 0.3535533905932738]
        @test Tₘ.P == [0.0 0.0; 0.0 0.0;;; 0.7071067811865476 0.0; 0.0 0.7071067811865476;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.7071067811865476 0.0; 0.0 0.7071067811865476;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.7071067811865476 0.0; 0.0 0.7071067811865476;;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0;;; 0.0 0.0; 0.0 0.0]
    end
end


#=
### *J₁-J₂ Chain* : *Power Method*
=#
bondD = 8
J1 = 1.0; J2 = 0.241167
Tₘ = model(J1J2Chain(), J1, J2)

βlist = [1.0, 2.0]
#Previous result with bond dimension 8, to_group=2
Op = zeros(2, 5)
Op[1, :] = [1.0, -0.7899920956097901, -0.18552679178607434, 0.6044653038237158]
Op[2, :] = [2.0, -0.5187573807353516, -0.29632701499366476, 0.4448607314833737] 
    
    
@testset "J1J2Chain_power_group" begin
    result_folder = "./J1J2Chain_power_group"
    to_group = 2
    power_evaluate(Tₘ, bondD, βlist; processor, result_folder, obsv_functions, to_group)
    Og = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Og ≈ Op  rtol = 1.e-5
end


@testset "J1J2Chain_power_shift" begin
    result_folder = "./J1J2Chain_power_shift"
    to_shift = 1.e-3
    power_evaluate(Tₘ, bondD, βlist; processor, result_folder, obsv_functions, to_shift)
    Os = readdlm(result_folder*"/obsvs.txt", skipstart = 1)
    @test Os ≈ Op  rtol = 1.e-5
end

