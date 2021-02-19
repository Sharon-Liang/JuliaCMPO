include("setup_struct.jl")

"""
cMPO for models
    (define functions in the future)
"""
X = [0. 1.; 1. 0]
Y = [0. -1im; 1im 0.]
Z = [1. 0.; 0. -1.]

"""
NN Transvers field Ising model
    H = ∑ J Zi Zj + ∑ Γ Xi
"""
J = 1. ; Γ = 1.
W = cMPO(Γ*X, √J*Z, √J*Z, zeros(2,2))
