using Arpack, LinearAlgebra
using Printf
using SparseArrays
using DelimitedFiles
using JLD, HDF5

function energy_density(J::Real, Γ::Real, N::Int, m::Int)
    k = m * π /N
    if k == 0
        ϵk = 2*(Γ - J)
    else
        ϵk = 2 * √(Γ^2 - 2*Γ*J*cos(k) + J^2)
    end
    return ϵk
end

struct pair
    subset::Vector{Int}
    Ptotal::Float64
end

function subsets(k::Int, Pmax::Number, Ps::Vector)
    """
        Given a set of several items with "prices" specified by the list Ps, find
        all subsets of length k whose total price do not exceed Pmax.
        Ps: energy of each k (sorted)
        subset : index of Ps
        Ptotal : Total energy of subset
    """
    Nelements = length(Ps)
    subset = zeros(Int, k); Ptotal = 0.
    result = [pair(subset, Ptotal)]
    for i = 1: k
        result_new = pair[]
        for p in result
            sum(p.subset) == 0 ? next_idx = 1 : next_idx = p.subset[i-1] + 1
            while next_idx + k - i < Nelements
                subset = copy(p.subset) 
                if sum(Ps[next_idx: next_idx+k-i]) <= Pmax - p.Ptotal
                    subset[i] = next_idx ;  Ptotal = p.Ptotal + Ps[next_idx]
                    append!(result_new, [pair(subset,Ptotal)])
                end
                next_idx += 1
            end
        end
        result = copy(result_new)
    end
    sort!(result, by = x -> x.Ptotal)
    indices, Ptotals = [x.subset for x in result], [x.Ptotal for x in result]
    return indices, Ptotals
end


function energy_level(J::Real, Γ::Real, N::Integer; 
    num::Int = 200, emax_diff::Real = 4, Nexcite_max::Int = 10)
    # N: Number of lattice sites
    # num: number of energy levels, default N
    # emax: maximum energy difference
    if N%4 != 0 @error "N must be a multiply of 4" end
    log2(num) > N ? num = 2^N : 
    Nexcite_max = max(N, Nexcite_max)

    # allowed momentum of even/odd sector
    m_even = [i for i in range(-(N-1), N-1, step = 2)]
    Em_even = [energy_density(J,Γ,N,m) for m in m_even] |> sort!
    m_odd = [i for i in range(-N+2, N, step = 2)]
    Em_odd = [energy_density(J,Γ,N,m) for m in m_odd] |> sort!
    
    # ground state of even sector
    gs_even = -0.5 * sum(Em_even)
    energy = [gs_even]

    # estimate k-limit
    Nresult = 1
    for i = 2:2:Nexcite_max
        indices, excite_even = subsets(i, emax_diff, Em_even)
        excite_even .+= gs_even
        length(indices) == 0 ? break : append!(energy, excite_even)
        println("lenth  of energy array= ", length(energy))
    end

    # odd sector
    gs_odd = -0.5 * sum(Em_odd) # not gs energy
    for i = 1:2:Nexcite_max # gs included
        indices, excite_odd = subsets(i, emax_diff, Em_odd)
        excite_odd .+= gs_odd
        length(indices) == 0 ? break : append!(energy, excite_odd)
        println("lenth  of energy array= ", length(energy))
    end

    sort!(energy)
    println("final length of energy arrag = ", length(energy))
    num = min(num, length(energy))
    println("final number of energy levels = ", num)
    return energy[1:num]    
end

"""check gamma/J """
J = [i for i in range(0.,1.,step = 0.05)]
sitenum = [16, 100]

for len in sitenum
    path1 = @sprintf "./data/exact/analytic-fixg%i.txt" len
    println(path1)
    println("Fix Γ, change J, sitenum = ", len)
    open(path1, "w") do file
        for j in J
            e = energy_level(j, 1., len, emax_diff = 8, Nexcite_max=len)
            writedlm(file,[j e'])
        end
    end
end