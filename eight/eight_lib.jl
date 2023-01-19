module eight_lib
export farthest_first_traversal

function farthest_first_traversal(points_arr::Array{<:Number}, num_centers::Int)::Array{<:Number}
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
            #curr_dist = minimum(sum((centers .- points_arr[i:i, :]) .^ 2, dims=2))
            curr_dist = minimum(sum((@. (centers - points_arr[i:i, :]) ^ 2), dims=2))
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

end
