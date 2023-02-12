module nine

export construct_trie
export prefix_trie_match
export lcs

function construct_trie(reads_vector)
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


function lcs(s1::AbstractString, s2::AbstractString)::String
    l, r, sub_len = 1, 0, 0
    for i in eachindex(s1)
        for j in i:length(s1)
            contains(s2, SubString(s1, i, j)) || break
            if sub_len ≤ j - i
                l, r = i, j
                sub_len = j - i
            end
        end
    end 
    return s1[l:r] 
end

end
