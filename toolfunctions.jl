# Tool fuctions
function symmetrize(A::AbstractArray)
    return (A + A')/2
end
