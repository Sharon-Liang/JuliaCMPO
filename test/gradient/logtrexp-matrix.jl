#logtrexp
using Test, JuliaCMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()
using FiniteDifferences, Zygote

ngradient(f, x) = grad(central_fdm(5, 1), f, x)
zgradient = Zygote.gradient

Nl, Nr = 3, 4

ψl, ψr = init_cmps(Nl), init_cmps(Nr)
Cmatrix = ψl * ψr
M = Cmatrix |> Matrix

t = rand()
∂t_M = x->logtrexp(x, M)
∂t_Cmatrix = x->logtrexp(x, Cmatrix)
∂t_Cmatrix_fed = x->logtrexp(x, Cmatrix, Full_ED)

∂M = x->logtrexp(t, x)
∂M_fed = x->logtrexp(t, x, Full_ED)

gradt_M_ndiff = ngradient(∂t_M, t)[1] 
gradt_Cmatrix_ndiff = ngradient(∂t_Cmatrix, t)[1] 
gradt_Cmatrix_fed_ndiff = ngradient(∂t_Cmatrix_fed, t)[1] 

@test gradt_Cmatrix_fed_ndiff ≈ gradt_Cmatrix_ndiff

gradM_ndiff = ngradient(∂M, M)[1] 
gradCmatrix_ndiff = ngradient(∂M, Cmatrix)[1] 
#gradCmatrix_fed_ndiff = ngradient(∂M_fed, Cmatrix)[1] 
gradCmatrix_fed_ndiff = grad(central_fdm(2, 1), ∂M_fed, Cmatrix)[1] 

@test ≈(gradCmatrix_fed_ndiff.ψl, gradCmatrix_ndiff.ψl; rtol = 1e-3, atol = 1e-5)
@test ≈(gradCmatrix_fed_ndiff.ψr, gradCmatrix_ndiff.ψr; rtol = 1e-5, atol = 1e-5)


gradt_M_zygote = zgradient(∂t_M, t)[1] 
gradt_Cmatrix_zygote = zgradient(∂t_Cmatrix, t)[1]
gradt_Cmatrix_fed_zygote = zgradient(∂t_Cmatrix_fed, t)[1]
gradM_zygote = zgradient(∂M, M)[1]
gradCmatrix_zygote = zgradient(∂M, Cmatrix)[1]
gradCmatrix_fed_zygote = zgradient(∂M_fed, Cmatrix)[1] 


@test ≈(gradt_M_zygote, gradt_M_ndiff; rtol = 1e-5, atol = 1e-5)
@test ≈(gradt_Cmatrix_zygote, gradt_Cmatrix_ndiff; rtol = 1e-5, atol = 1e-5)
    
@test ≈(gradM_zygote, gradM_ndiff; rtol = 1e-5, atol = 1e-5)

@test ≈(gradCmatrix_zygote.ψl, gradCmatrix_ndiff.ψl; rtol = 1e-5, atol = 1e-5)
@test ≈(gradCmatrix_zygote.ψr, gradCmatrix_ndiff.ψr; rtol = 1e-5, atol = 1e-5)