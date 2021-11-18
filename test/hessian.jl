using cMPO, Test
using Random; Random.seed!()
using Zygote, ForwardDiff, FiniteDiff

function hessiancheck(f, xs)
    h_zygote = Zygote.hessian(f, xs)
    h_finite_difference = FiniteDiff.finite_difference_hessian(f, xs)
    return all(isapprox.(h_zygote, h_finite_difference; rtol = 1e-5, atol = 1e-5))
end
hessiantest(f, xs::AbstractVector) = hessiancheck(xs -> sum(sin.(f(xs))), xs)

function hessiancheck_fdiff(f, xs)
    h_zygote = ForwardDiff.hessian(f, xs)
    h_finite_difference = FiniteDiff.finite_difference_hessian(f, xs)
    return all(isapprox.(h_zygote, h_finite_difference; rtol = 1e-5, atol = 1e-5))
end
hessiantest_fdiff(f, xs::AbstractArray) = hessiancheck_fdiff(xs -> sum(sin.(f(xs))), xs)


function test_prod(sl::AbstractVector, dl::Tuple, sr::AbstractVector, dr::Tuple)
    sl = tocmps(sl, dl)
    sr = tocmps(sr, dr)
    sl * sr
end

function test_prod(sl::AbstractVector, dl::Tuple, o::CMPO; 
    sr::AbstractVector, dr::Tuple)
    sl = tocmps(sl, dl)
    sr = tocmps(sr, dr)
    sl * o * sr
end

function test_prod(o::CMPO, sr::AbstractVector, dr::Tuple; 
    sl::AbstractVector, dl::Tuple)
    sr = tocmps(sr, dr)
    sl = tocmps(sl, dl)
    sl * (o * sr)  #vs: sl * o * sr 先计算哪个看起来计算图有区别,左右不同时会有区别
end

@testset "*: zygote" begin
    χ = 2
    sl, dl = tovector(init_cmps(χ))
    sr, dr = tovector(init_cmps(χ))
    g = rand(1)[1]
    o = TFIsing(1.0, g)

    @test hessiantest(x -> test_prod(x, dl, sr, dr) , sl)
    @test hessiantest(x -> test_prod(sl, dl, x, dr) , sr)
    @test hessiantest(x -> test_prod(x,dl,o; sr=sr, dr=dr), sl)
    @test hessiantest(x -> test_prod(o,x,dr; sl=sl, dl=dl), sr)

    @test hessiantest(x -> test_prod(x, dl, x, dl) , sl)
    @test hessiantest(x -> test_prod(x,dl,o;sr=x, dr=dl), sl)
    @test hessiantest(x -> test_prod(o,x,dr;sl=x, dl=dr), sr)
end

@testset "*: forwarddiff" begin
    χ = 2
    sl, dl = tovector(init_cmps(χ))
    sr, dr = tovector(init_cmps(χ))
    g = rand(1)[1]
    o = TFIsing(1.0, g)

    @test hessiantest_fdiff(x -> test_prod(x, dl, sr, dr) , sl)
    @test hessiantest_fdiff(x -> test_prod(sl, dl, x, dr) , sr)
    @test hessiantest_fdiff(x -> test_prod(x,dl,o; sr=sr, dr=dr), sl)
    @test hessiantest_fdiff(x -> test_prod(o,x,dr; sl=sl, dl=dl), sr)
    @test hessiantest_fdiff(x -> test_prod(x,dl,o,sr,dr), sl)

    @test hessiantest_fdiff(x -> test_prod(x, dl, x, dl) , sl)
    @test hessiantest_fdiff(x -> test_prod(x,dl,o;sr=x, dr=dl), sl)
    @test hessiantest_fdiff(x -> test_prod(o,x,dr;sl=x, dl=dr), sr)
end

@testset "logtrexp" begin
    D = 2
    A = rand(D, D) |> symmetrize
    @test hessiantest(x->logtrexp(reshape(x,D,D)), vec(A))  
    @test hessiantest_fdiff(x->logtrexp(reshape(x,D,D)), vec(A))  
end

@testset "free energy: TFIsing model" begin
    β = 20.0
    χ = 2
    g = rand(1)[1]
    w = TFIsing(1.,g)
    ψ, dim = init_cmps(χ) |> tovector
    Zygote.hessian(x->free_energy(x, dim, w, β), ψ)
    @test hessiancheck(x->free_energy(x, dim, w, β), ψ)
    @test hessiancheck_fdiff(x->free_energy(x, dim, w, β), ψ)
end