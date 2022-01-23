using Arpack, LinearAlgebra, SparseArrays
using DelimitedFiles, JLD, HDF5
using Printf
using Yao

function TFIsing(J::Real, Γ::Real, N::Integer; PBC::Bool=true)
    J = sqrt(J)
    hloc = sum([repeat(N, Γ*X, i) for i=1:N]) #Transverse field terms
    hint = sum([repeat(N, J*Z, (i, i+1)) for i = 1:N-1])
    if PBC 
        hint += repeat(N, J*Z, (N,1))
    end
    return hloc + hint
end

sitenum = [7,8,10,16]
eignum = 128


J = [i for  i in range(0,1, step=0.1)]
gamma = [i for i in range(0,1,step=0.1)]
for len in sitenum
    path1 = @sprintf "../data/exact/exact-fixg%i.txt" len
    path2 = @sprintf "../data/exact/exact-fixj%i.txt" len

    println("Fix Γ, change J, sitenum = ", len)
    open(path1, "w") do file
        for j in J
            h = TFIsing(j,1.0,len)
            mh = Matrix(h) |> real |> Hermitian |> sparse
            if len == 7
                e = eigvals(mh)
            else
                e = eigs(mh, nev=eignum, which = :SR)[1]
            end
            e = e ./ len
            writedlm(file,[j e'])
        end
    end


    println("Fix J, change Γ, sitenum = ", len)
    open(path2, "w") do file
        for Γ in gamma
            h = TFIsing(1.0,Γ,len)
            mh = Matrix(h) |> real |> Hermitian |> sparse
            if len == 7
                e = eigvals(mh)
            else
                e = eigs(mh, nev=eignum, which = :SR)[1]
            end
            writedlm(file,[Γ e'])
        end
    end
end
            