#module Models

#using LinearAlgebra

#include("ToolFunctions.jl")
#include("SetupStruct.jl")

"""
cMPO for models
    (define functions in the future)
"""

"""
Pauli matrices
"""
X = [0. 1.; 1. 0]
Y = [0. -1im; 1im 0.]
Z = [1. 0.; 0. -1.]


"""
NN Transvers field Ising model
    H = ∑ J Zi Zj + ∑ Γ Xi
"""
J = 1. ; Γ = 1.
W = cmpo(Γ*X, √J*Z, √J*Z, zeros(2,2))

#end  # modul Models
