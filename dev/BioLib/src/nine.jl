module nine

import DataStructures: DefaultDict

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
        Return index of beginning of text pattern from original trie text. Return empty string if trie don't contain text string
    =#
    symb = first(text)
    v = 0 #current node in trie. Assume that root node have label 0
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
    Find longest common substring between two strings using dynamic programming
    by checking lengths of all common prefixes
    Input data:
    -----------
        str1 - first string
        str2 - second string
    Output data:
    -----------
        longest common substring between str1 and str2
    =#
    lcs_arr = zeros(Int32, map(x -> length(x) + 1, (str1, str2)))
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


"""
Naive implementation of burrowhs-wheels transform
Input string should not contain stx and edx symbols
stx - unicode start of text symbol
edx - unicode end of text symbol
"""
function bwt_naive(inp_str::AbstractString; stx="\x02", edx="\x03") :: AbstractString
    #if match(r"\x02|\x03", inp_str) !== nothing
    #    throw("String for Burrows-Wheeler input cannot contain STX or ETX")
    #end
    inp_str = stx * inp_str * edx
    return String([t[end] for t in sort([circshift([c for c in inp_str], n) for n in 0:length(inp_str)-1])])
end


function ibwt_naive(bwt_str::AbstractString)
    char_vector_last_column = Vector{Char}(bwt_str)
    char_vector_first_column = String(sort(char_vector_last_column))
    skip_count = 0#DefaultDict{Char, Int}(0)
    curr_char = '\$'
    res = Vector{Char}()
    for _ in 1:length(bwt_str)
        idx = findfirst(curr_char, bwt_str)
        while skip_count != 0
            for i in bwt_str[idx + 1:end]
                idx += 1
                if  curr_char == i
                    break
                end
            end
            skip_count -= 1
        end

        curr_char = char_vector_first_column[idx]
        push!(res, curr_char)

        skip_count = 0
        skip_char_idx = idx - 1
        while skip_char_idx != 0 && char_vector_first_column[skip_char_idx] == curr_char
            skip_count += 1
            skip_char_idx -= 1
        end
    end
    return String(res)
end

function compute_last2first(bwt_str::AbstractString)
    skip_count = DefaultDict{Char, Int}(0)
    bw_first = sort(Vector{Char}(bwt_str))
    res = zeros(Int, length(bw_first))
    for (idx, c) in enumerate(bwt_str)
        res[idx] = searchsortedfirst(bw_first, c) + skip_count[c]
        skip_count[c] += 1
    end
    return res
end

function bwt_matching(bwt_last::AbstractString, pattern::Vector{Char}, last2first::Vector{Int})
    top = 1
    bottom = length(bwt_last)
    while top <= bottom
        if length(pattern) > 0
            curr_symb = pop!(pattern)

            top_idx = - 1
            bottom_idx = -1
            top_not_found_flag = true
            for (idx, i) in enumerate(bwt_last[top:bottom])
                if i == curr_symb 
                    if top_not_found_flag
                        top_idx = idx + top - 1
                        bottom_idx = idx + top - 1
                        top_not_found_flag = false
                    else
                        bottom_idx = idx + top - 1
                    end
                end
            end

            if top_idx == -1
                return 0
            end
            top = last2first[top_idx]
            bottom = last2first[bottom_idx]
        else
            return bottom - top + 1
        end
    end
    return -1
end


end
