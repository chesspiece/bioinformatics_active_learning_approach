module nine

export construct_trie
export prefix_trie_match
export lcs

function construct_trie(reads_vector::Vector{<:AbstractString})::Dict{Int32,Dict{Char,Int32}}
    #=
    Construct trie data strcutures from vector of strings
    Input:
    -----
        reads_vector - vector which have srings of dna as elements
    Outputs:
    --------
        Trie data strcuture which was build using strings from reads_vector
    =#
    trie = Dict{Int32,Dict{Char,Int32}}(0 => Dict())
    num_vertexes = 1
    for read in reads_vector
        curr_node = 0
        for nucleotide in read
            if nucleotide in keys(trie[curr_node])
                curr_node = trie[curr_node][nucleotide]
            else
                trie[curr_node][nucleotide] = num_vertexes
                trie[num_vertexes] = Dict()
                curr_node = num_vertexes
                num_vertexes += 1
            end
        end
    end

    return trie
end

function prefix_trie_match(text::String, trie::Dict{Int32,Dict{Char,Int32}})::Tuple{Bool,AbstractString}
    #=
    Check if substring is presented in the trie
    Input data:
    ----------
        text - check if this string is contaned in trie
        trie - trie data strcuture which represent vector of string in good for searcing form
    Output data:
    -----------
        Return true if text is contained in trie, False otherwise.
        Return index of starteing of text pattern from original trie text. Return empty string if trie don't contain text string
    =#
    symb = first(text)
    v = 0
    tl = 1
    while true
        if length(trie[v]) == 0
            return true, text[1:tl-1]
        elseif symb in keys(trie[v])
            symb = text[tl]
            v = trie[v][symb]
            tl += 1
            if tl <= length(text)
                symb = text[tl]
            elseif length(trie[v]) == 0
                return true, text[1:tl-1]
            else
                return false, ""
            end
        else
            return false, ""
        end
    end
end


function lcs(str1::AbstractString, str2::AbstractString)::AbstractString
    #=
    Find longest common substring between two strings using dynamic programming by checking lengths of all common prefixes
    =#
    lcs_arr = zeros(Int32, map(x -> x + 1, (length(str1), length(str2))))
    #lcs_arr[1:1, :] .= 0
    #lcs_arr[:, 1:1] .= 0
    pos = 0
    mx_ln = 0

    for i in 1:length(str1)
        for j in 1:length(str2)
            if str1[i] == str2[j]
                lcs_arr[i + 1, j + 1] = lcs_arr[i, j] + 1
                if mx_ln < lcs_arr[i + 1, j + 1] 
                    mx_ln = lcs_arr[i + 1, j + 1]
                    pos = j
                end
            end
        end
    end
    return str2[pos-mx_ln+1:pos]
end


end
