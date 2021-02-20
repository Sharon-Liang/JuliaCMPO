module ToolFunctions

using LinearAlgebra

function symmetrize(A::AbstractArray)
    return (A + A')/2
end

function tr_exp(A::AbstractMatrix, β::Real)
    vals = eigvals(A)
    res = 0.
    for i = 1:length(vals)
        res += exp(-β*vals[i])
    end
    res = real(res)
    if res < 0.
        println("Warnning: Negative tr_exp!")
    end
    return res
end


end  # module ToolFunctions
