using Flux
using HCubature   

const x_max = 1
actual(x) = x
function actual_int(k)
    res, _ = hquadrature(x -> actual(x)* k, 0, x_max)
    return res
end

k_train, k_test = hcat(1:5...), hcat(6:10...)
y_train, y_test = actual_int.(k_train), actual_int.(k_test)

predict = Dense(1,1)
parameters = params(predict)
fpred(x) = predict([x])[1]

function int_func(k)
    N = 1000
    dx = x_max / N
    xrange = hcat([i for i=0:dx:x_max]...)
    res = sum(fpred.(xrange) .* k) * dx
    return res
end


function loss(k, y)
    ypred = int_func.(k)
    return Flux.mse(ypred, y)
end

l1 = loss(k_train, y_train)

gradient(()->loss(k_train, y_train), params(predict))


using Flux: train!
opt = Descent(0.05)
data = [(k_train, y_train)]
train!(loss, parameters, data, opt)

for epoch in 1:200
    train!(loss, parameters, data, opt)
end