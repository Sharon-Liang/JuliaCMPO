using Test
using Random ; Random.seed!()
using Zygote
using Optim
using cMPO
@test true

"""
Gradient test
1. Zygote gradient test
2. Optim autodiff gradient test
"""
