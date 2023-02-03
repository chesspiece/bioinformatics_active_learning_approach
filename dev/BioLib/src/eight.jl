module eight
export farthest_first_traversal
export compute_mse_distance
export lloyd_alg

function farthest_first_traversal(points_arr::Array{<:Number,2}, num_centers::Int)::Array{<:Number,2}
    _, d2 = size(points_arr)
    centers = zeros((num_centers, d2))
    centers[1, :] = points_arr[1, :] # let's initialize first center
    save_idx = [1]
    for centers_row in 2:num_centers
        elm = nothing
        dist = -Inf
        for i in axes(points_arr, 1)
            if i in save_idx
                continue
            end
            curr_dist = minimum(sum((@. (centers - points_arr[i:i, :])^2),
                dims=2))
            if curr_dist > dist
                dist = curr_dist
                elm = i
            end
        end
        centers[centers_row, :] = points_arr[elm, :]
        push!(save_idx, elm)
    end
    return centers
end


function compute_mse_distance(centers::Array{<:Number,2}, datapoints::Array{<:Number,2})::Number
    num_points, _ = size(datapoints)
    sum_dist = 0
    for datapoint in axes(datapoints, 1)
        sum_dist += minimum(sum((@. (centers - datapoints[datapoint:datapoint, :])^2),
            dims=2))
    end
    return sum_dist / num_points
end


function lloyd_alg(points_arr::Array{<:Number,2}, num_centers::Int)::Array{<:Number,2}
    _, d2 = size(points_arr)
    centers_idx = [i for i in 1:num_centers]
    centers = points_arr[centers_idx, :]
    while true
        new_centers = zeros((num_centers, d2))
        new_centers_count = zeros(num_centers)
        for i in axes(points_arr, 1)
            curr_dist = sum((@. (centers - points_arr[i:i, :])^2), dims=2)
            curr_dist = dropdims(curr_dist, dims=tuple(findall(size(curr_dist) .== 1)...))
            idx = argmin(curr_dist)
            new_centers[idx:idx, :] = new_centers[idx:idx, :] .+ points_arr[i:i, :]
            new_centers_count[idx] = new_centers_count[idx] + 1
        end
        new_centers = new_centers ./ new_centers_count
        diff = maximum(abs.(new_centers - centers))
        centers = new_centers
        if diff < 1e-14
            break
        end
    end
    return centers
end

end

