module eight_runners

using BioLib.eight
using Printf

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


function second_task()
    k, m, centers_arr, points_arr = open("./input/input_8_2.txt", "r") do f
        k, m = parse.(Int, split(readline(f), " "))
        centers_arr = Vector{Array{Float64}}(undef, 0)
        points_arr = Vector{Array{Float64}}(undef, 0)
        for _ in 1:k
            push!(centers_arr, parse.(Float64, split(readline(f), " ")))
        end
        readline(f)
        for coord in eachline(f)
            push!(points_arr, parse.(Float64, split(coord, " ")))
        end
        return k, m, centers_arr, points_arr
    end

    res = compute_mse_distance(mapreduce(permutedims, vcat, centers_arr),
        mapreduce(permutedims, vcat, points_arr))
    open("output.txt", "w") do f
        println(f, "$(@sprintf("%.3f", res))")
    end
    # return points_arr, centers_arr, res
end

function third_task()
    k, m, points_arr = open("./input/input_8_3.txt", "r") do f
        k, m = parse.(Int, split(readline(f), " "))
        points_arr = Vector{Array{Float64}}(undef, 0)
        for coord in eachline(f)
            push!(points_arr, parse.(Float64, split(coord, " ")))
        end
        return k, m, points_arr
    end

    res = lloyd_alg(mapreduce(permutedims, vcat, points_arr), k)
    open("output.txt", "w") do f
        for i in axes(res, 1)
            s = join(map(x -> @sprintf("%.3f", x), res[i, :]), " ")
            println(f, s)
        end
    end
    # return points_arr, centers_arr, res
end


end
