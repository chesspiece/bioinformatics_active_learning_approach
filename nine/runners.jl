module nine_runners

import BioLib.nine: construct_trie, prefix_trie_match, lcs, bwt_naive, ibwt_naive,
       compute_last2first, bwt_matching
import SfTree: build_sf_tree, get_edges_names, longest_repeated_substring
using SuffixArrays
using Printf


function first_task()
    reads_arr = open("./input/input_9_1.txt", "r") do f
        reads_arr = split(readline(f), " ")
        return reads_arr
    end

    trie = construct_trie(reads_arr)
    needed_vertexes = filter(n -> length(trie[n]) > 0, keys(trie))

    open("output.txt", "w") do f
        for v in needed_vertexes
            for (edge, s_v) in trie[v]
                s_v = string(s_v)
                println(f, join([v, s_v, edge], " "))
            end
        end
    end
    #return reads_arr, construct_trie(reads_arr)
end


function second_task()
    reference_dna, reads_arr = open("./input/input_9_2.txt", "r") do f
        reference_dna = readline(f)
        reads_arr = split(readline(f), " ")
        return reference_dna, reads_arr
    end

    trie = construct_trie(reads_arr)
    res::Dict{String,Vector{Int}} = Dict()
    for i in reads_arr
        res[i] = []
    end

    for i in eachindex(reference_dna)
        flag, str_key = prefix_trie_match(reference_dna[i:end], trie)
        if flag
            push!(res[str_key], i - 1)
        end
    end

    open("output.txt", "w") do f
        for key in keys(res)
            println(f, key * ": " * join(string.(res[key]), " "))
        end
    end
end


function third_task()
    dna = open("./input/input_9_3.txt", "r") do f
        dna = readline(f)
        return dna
    end
    sf_tree = build_sf_tree(dna)
    res = get_edges_names(sf_tree)
    open("output.txt", "w") do f
        println(f, join(res, " "))
    end
end


function fourth_task()
    dna = open("./input/input_9_4.txt", "r") do f
        dna = readline(f)
        return dna
    end
    sf_tree = build_sf_tree(dna)
    res = get_edges_names(sf_tree)
    open("output.txt", "w") do f
        println(f, join(res, " "))
    end
end


function fifth_task()
    dna = open("./input/input_9_5.txt", "r") do f
        dna = readline(f)
        dna *= "#"
        return dna
    end
    sf_tree = build_sf_tree(dna)
    open("output.txt", "w") do f
        println(f, longest_repeated_substring(sf_tree))
    end
end


function sixth_task()
    dna1, dna2 = open("./input/input_9_6.txt", "r") do f
        dna1 = readline(f)
        dna2 = readline(f)
        return dna1, dna2
    end

    open("output.txt", "w") do f
        println(f, lcs(dna1, dna2))
    end
end


function seventh_task()
    dna = open("./input/input_9_6.txt", "r") do f
        dna1 = readline(f)
        dna2 = readline(f)
        dna = dna1 * "#" * dna2 * "!"
        return dna
    end
    sf_tree = build_sf_tree(dna)
    open("output.txt", "w") do f
        println(f, longest_repeated_substring(sf_tree))
    end
end


function eighth_task()
    dna = open("./input/input_9_8.txt", "r") do f
        dna = readline(f)
        return dna
    end

    open("output.txt", "w") do f
        println(f, join(string.(map(x -> x - 1, suffixsort(dna)), base=10), " "))
    end
end


function nineth_task()
    dna = open("./input/input_9_9.txt", "r") do f
        dna = readline(f)
        return dna
    end
    open("output.txt", "w") do f
        println(f, bwt_naive(dna, stx="", edx=""))
    end
end


function tenth_task()
    dna = open("./input/input_9_10.txt", "r") do f
        dna = readline(f)
        return dna
    end
    open("output.txt", "w") do f
        println(f, ibwt_naive(dna))
    end
end


function eleventh_task()
    dna_bwt, patterns = open("./input/input_9_11.txt", "r") do f
        dna_bwt = readline(f)
        patterns = split(readline(f), " ")
        return dna_bwt, patterns
    end
    l2f = compute_last2first(dna_bwt)
    res = Vector{Int}()
    for pat in patterns
        push!(res, bwt_matching(dna_bwt, Vector{Char}(pat), l2f))
    end
    open("output.txt", "w") do f
        println(f, join(res, " "))
    end
end


end
