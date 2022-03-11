#module SaveAndLoad
function saveCMPS(path::String, ψ::CMPS, dict::Dict = Dict())
    h5open(path, "w") do file
        file["Q"] = ψ.Q
        file["R"] = ψ.R
        for key in keys(dict) 
            file[key] = dict[key] 
        end
    end
end

function readCMPS(path::String; python::Bool = false)
    h5open(path, "r") do file
        Q = read(file["Q"]) 
        R = read(file["R"])
        if python
            Q = Q |> diagm
            R = ein"ijk->jik"(R)
        end
        return CMPS(Q, R)
    end
end

#end    # module SaveAndLoad