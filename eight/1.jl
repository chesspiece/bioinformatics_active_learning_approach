module eight_runners

include("./eight_lib.jl")
using .eight_lib

function first_task()
    k, m, points_arr = open("./input/input_8_1.txt", "r") do f
        k, m = parse.(Int, split(readline(f), " "))
        points_arr = Vector{Array{Float64}}(undef, 0)
        for coord in eachline(f)
            push!(points_arr, parse.(Float64, split(coord, " ")))
        end
        return k, m, points_arr
    end

    res = farthest_first_traversal(mapreduce(permutedims, vcat, points_arr), k)
    open("output.txt", "w") do f
       for i in axes(res, 1)
           s = join(string.(res[i, :]), " ")
           println(f, s)
       end
    end
end

end
