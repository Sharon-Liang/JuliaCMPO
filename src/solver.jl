
TensorUnion = Union{AbstractCTensor, CMPSMatrix}

"""
    `cpu_solver` for CTensors and CMPSMatrix
"""
function cpu_solver(f::Function, M::T...) where T<:TensorUnion
    #M = Tuple(CTensor(x) for x in M)
    M = map(x->CTensor(x), M)
    return f(M...)
end

"""
    `gpu_solver` for CTensors and CMPSMatrix
"""
function gpu_solver(f::Function, M::T...) where T<:TensorUnion
    #M = Tuple(CuCTensor(x) for x in M)
    M = map(x->CuCTensor(x), M)
    return f(M...)
end
