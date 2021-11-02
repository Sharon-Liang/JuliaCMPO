using Pkg
Pkg.activate("../")
using cMPO
using Optim
using Printf

function converttime(t::Number)
    (m, s) = divrem(t, 60)
    (h, m) = divrem(m, 60)
    (d, h) = divrem(m, 24)
    if d == 0
        if h == 0
            if m == 0 return @sprintf "time=%.3fs" s
            else return @sprintf "time=%imin %.3f s" m s
            end
        else
            return @sprintf "time=%ih %imin %.3fs" h m s
        end
    else
        return @sprintf "time=%id %ih %imin %.3fs" d h m s
    end
end

println("2021-11-02:test NewtonTrustRegion time")
println("χ = 8; β = 1.0")
χ = 8; β = 1.0
w = TFIsing(1.0, 1.0);
arr = init_cmps(χ,w) |> toarray;
(vec, dim) = init_cmps(χ,w) |> tovector;
f1 = x -> free_energy(x, w, β);
f2 = x -> free_energy(x, dim, w, β);
g1 = gradient_function(f1);

g2 = gradient_function(f2);
h = hessian_function(f2);


t1 = time();
op1 = optimize(f1, g1, arr, LBFGS(), Optim.Options(iterations = 10000));
t2 = time();
println("LBFGS ", converttime(t2-t1))


t1 = time();
op2 = optimize(f2, g2, h, vec, NewtonTrustRegion());
t2 = time();
println("NewtonTrustRegion ", converttime(t2-t1))
