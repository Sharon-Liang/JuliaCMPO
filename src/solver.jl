function cpu_solver(f::Function, M::AbstractArray...)
    return f(M...)
end

function gpu_solver(f::Function, M::AbstractArray...)
    M = Tuple(convert(CuArray, x) for x in M)
    eltype(M) <: CuArray ? (return f(M...)) : error("CuArray conversion fails")
end

function cpu_solver(f::Function, M::AbstractCTensor...)
    M = Tuple(CTensor(x) for x in M)
    eltype(M) <: Union{CMPS, CMPO} ? (return f(M...)) : error("CTensor conversion fails")
end

function gpu_solver(f::Function, M::AbstractCTensor...)
    M = Tuple(CuCTensor(x) for x in M)
    eltype(M) <: Union{CuCMPS, CuCMPO} ? (return f(M...)) : error("CuCTensor conversion fails")
end
