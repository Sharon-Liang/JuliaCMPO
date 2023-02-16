#=
### *Processor Supported*
=#
"""
    Processor CPU GPU

    Enumerate processors: GPU and GPU
"""
@enum Processor CPU GPU


#=
### *Global Constants* : *Literal Strings*
=#

"""
    __LIBNAME__

Name of this julia toolkit.

See also: [`__VERSION__`](@ref).
"""
const __LIBNAME__ = "JuliaCMPO"

#=
"""
    __VERSION__

Version of this julia toolkit.

See also: [`__RELEASE__`](@ref).
"""
const __VERSION__ = v"1.3.6-devel.230129"

"""
    __RELEASE__

Release date of this julia toolkit.

See also: [`__AUTHORS__`](@ref).
"""
const __RELEASE__ = "2023/01"
=#


#=
*Remarks* :

The elements of the Array `__AUTHORS__` should be a `NamedTuple` object,
such as:

```julia
(name = "author's name", email = "author's email")
```
=#

"""
    __AUTHORS__

Core authors of this julia toolkit.

See also: [`__LIBNAME__`](@ref).
"""
const __AUTHORS__ = [(name = "Shuang Liang", email = "sliang@iphy.ac.cn")]

"""
    authors()

Print authors / contributors of the `NevanlinnaAC` toolkit.

See also: [`__AUTHORS__`](@ref).
"""
function authors()
    println("Authors (Until $__RELEASE__):")
    for a in __AUTHORS__
        println("  $(a.name) (email: $(a.email))")
    end
end