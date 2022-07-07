function cpu_solver(f::Function, M::AbstractArray...)
    M = map(x->convert(Array, x), M)
    return f(M...)
end

function gpu_solver(f::Function, M::AbstractArray...)
    M = map(x->convert(CuArray, x), M)
    return f(M...)
end

TensorUnion = Union{AbstractCTensor, CMPSMatrix}
function cpu_solver(f::Function, M::T...) where T<:TensorUnion
    #M = Tuple(CTensor(x) for x in M)
    M = map(x->CTensor(x), M)
    return f(M...)
end

function gpu_solver(f::Function, M::T...) where T<:TensorUnion
    #M = Tuple(CuCTensor(x) for x in M)
    M = map(x->CuCTensor(x), M)
    return f(M...)
end
