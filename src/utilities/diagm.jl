#module utilities
"""
    diagm(v) function for CuVector
"""
LinearAlgebra.diagm(v::CuVector) = ein"i->ii"(v)
#end  # module utilities
