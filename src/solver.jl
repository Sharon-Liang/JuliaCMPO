
TensorUnion = Union{AbstractCTensor, CMPSMatrix}
"""
    cpu_solver(x::Union{AbstractCTensor, CMPSMatrix})
    cpu_solver(f::Function, x::Union{AbstractCTensor, CMPSMatrix})

Use CPU to solve ``f(x)``
"""
cpu_solver(x::TensorUnion) = CTensor(x)
cpu_solver(f::Function, x::TensorUnion) = f(CTensor(x))


"""
    gpu_solver(x::Union{AbstractCTensor, CMPSMatrix})
    gpu_solver(f::Function, x::Union{AbstractCTensor, CMPSMatrix})

Use GPU to solve ``f(x)``
"""
gpu_solver(x::TensorUnion) = CuCTensor(x)
gpu_solver(f::Function, x::TensorUnion) = f(CuCTensor(x))
