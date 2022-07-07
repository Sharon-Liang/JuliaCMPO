#logtrexp
using Test, JuliaCMPO
using LinearAlgebra, Parameters, FiniteTLanczos
using Random; Random.seed!()
using FiniteDifferences, Zygote

Nl, Nr = 3, 4

ψl, ψr = init_cmps(Nl), init_cmps(Nr)
Cmatrix = ψl * ψr
M = Cmatrix |> Matrix

t = rand()
∂t_M = x->logtrexp(x, M)
∂t_Cmatrix = x->logtrexp(x, Cmatrix)
∂M = x->logtrexp(t, x)

gradt_M_ndiff = ngradient(∂t_M, t)[1] 
gradt_Cmatrix_ndiff = ngradient(∂t_Cmatrix, t)[1] 

gradM_ndiff = ngradient(∂M, M)[1] 
gradCmatrix_ndiff = ngradient(∂M, Cmatrix)[1] 




gradt_M_zygote = zgradient(∂t_M, t)[1] 
gradt_Cmatrix_zygote = zgradient(∂t_Cmatrix, t)[1]
gradM_zygote = zgradient(∂M, M)[1]
gradCmatrix_zygote = zgradient(∂M, Cmatrix)[1]
    


@test ≈(gradt_M_zygote, gradt_M_ndiff; rtol = 1e-5, atol = 1e-5)
@test ≈(gradt_Cmatrix_zygote, gradt_Cmatrix_ndiff; rtol = 1e-5, atol = 1e-5)
    
@test ≈(gradM_zygote, gradM_ndiff; rtol = 1e-5, atol = 1e-5)

@test ≈(gradCmatrix_zygote.ψl, gradCmatrix_ndiff.ψl; rtol = 1e-5, atol = 1e-5)
@test ≈(gradCmatrix_zygote.ψr, gradCmatrix_ndiff.ψr; rtol = 1e-5, atol = 1e-5)