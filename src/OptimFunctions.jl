#module OptimFunctions
#Reference: https://github.com/baggepinnen/FluxOptTools.jl/blob/master/src/FluxOptTools.jl
veclength(grads::Zygote.Grads) = sum(length(grads[p]) for p in grads.params)
veclength(params::Zygote.Params) = sum(length, params.params)

Base.zeros(grads::Zygote.Grads) = zeros(veclength(grads))
Base.zeros(pars::Zygote.Params) = zeros(veclength(pars))

"""
    optim_function(loss, pars): Return two functions(loss function, gradient function)
    and p0, a vectorized version of pars.
"""
function optim_functions(loss, pars::Zygote.Params)
    grads = Zygote.gradient(loss, pars)
    p0 = zeros(pars)
    copy!(p0, pars)

    gradient_function = function (g,w)
        copy!(pars, w)
        grads = Zygote.gradient(loss, pars)
        copy!(g, grads)
    end

    loss_function = function (w)
        copy!(pars, w)
        loss()
    end
    return p0, loss_function, gradient_function
end


#end #module OptimFunctions





"""
    Generate gradient function which is used in Optim.jl

function gradfunc_gen(f::Function)
    function func(val::Array{T}, var::Array{T}) where T<:Number
        val[1:end] = gradient(f,var)[1][1:end]
    end
    return func
end

function hessfunc_gen(f::Function)
    function func(val::Matrix{T}, var::Vector{T}) where T<:Number
        val[1:end] = hessian(f,var)[1:end]
    end
    return func
end
"""