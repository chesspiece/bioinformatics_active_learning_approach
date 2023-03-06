module nine_runners

import BioLib.nine: construct_trie, prefix_trie_match, lcs
using SuffixTrees
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


function fourth_task()
    dna = open("./input/input_9_4.txt", "r") do f
        dna = readline(f)
        return dna
    end

    open("output.txt", "w") do f
        println(f, getlongestrepeatedsubstring(SuffixTree(dna)))
    end
end


function fifth_task()
    dna1, dna2 = open("./input/input_9_5.txt", "r") do f
        dna1 = readline(f)
        dna2 = readline(f)
        return dna1, dna2
    end

    open("output.txt", "w") do f
        println(f, lcs(dna1, dna2))
    end
end


function seventh_task()
    dna = open("./input/input_9_7.txt", "r") do f
        dna = readline(f)
        return dna
    end

    open("output.txt", "w") do f
        println(f, join(string.(map(x -> x - 1, suffixsort(dna)), base=10), " "))
    end
end

end
