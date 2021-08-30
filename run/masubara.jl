using cMPO
using LinearAlgebra
using DelimitedFiles
using JLD, HDF5
using Printf

g = 1.0
path = @sprintf "./data/g_%.1f.jld" g
d1 = load(path);
β = 10.0
key = string(β)

w = TFIsing(1.,g)
ψ = d1[key][2] |> tocmps
z = make_operator(pauli(:z), ψ)
Nfreq = 20

path1 = @sprintf "./data/Acmpo.txt"
path2 = @sprintf "./data/input.txt" # frequency real_part imag_part \n, increasing masubara frenquancies

open(path1, "w") do file
end

open(path2, "w") do file
end

for i in range(-10,10,step=0.05)
    open(path1, "a") do file1
        writedlm(file1,[i spectral_function(i,z,z,ψ,w,β,η=0.001)])
    end
end

for i = 1:Nfreq
    freq = Masubara_freq(i,β) |> imag
    gn = Masubara_GF(i,z,z,ψ,w,β)
    open(path2, "a") do file2
        writedlm(file2,[freq real(gn) imag(gn)])
    end
end
