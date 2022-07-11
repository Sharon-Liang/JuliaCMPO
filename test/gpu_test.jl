using JuliaCMPO, Test
using CUDA; CUDA.allowscalar(false)
solver = gpu_solver