module ten_runners

using BioLib.ten
using Printf


function first_task()
    read, states, prob_dict = open("./input/input_10_1.txt", "r") do f
        read = readline(f)
        readline(f)
        states = split(strip(readline(f)))
        readline(f)
        cols = split(strip(replace(readline(f), r" +" => " ")))
        #prob_dict = Dict{SubString{String}, Dict{SubString{String}, Float64}}()
        prob_dict = Dict()
        for cls in cols
            prob_dict[cls] = Dict()
        end
        for line in eachline(f)
            str_lin = split(strip(replace(line, r" +" => " ")))
            cur_c = str_lin[1]
            str_lin = parse.(Float64, str_lin[2:end]) 
            for (cls, prb) in zip(cols, str_lin)
                prob_dict[cur_c][cls] = prb
            end
        end
        return read, states, prob_dict
    end

    res = 0.5
    c_l = read[1]
    for c in read[2:end]
        res *= prob_dict[string(c_l)][string(c)]
        c_l = c
    end

    open("output.txt", "w") do f
        println(f, string(res))
    end
end


function second_task()
    str_c, read, states, prob_dict = open("./input/input_10_2.txt", "r") do f
        str_c = readline(f)
        readline(f)
        alphabet = split(strip(readline(f)))
        readline(f)
        read = readline(f)
        readline(f)
        states = split(strip(readline(f)))
        readline(f)

        cols = split(strip(replace(readline(f), r" +" => " ")))
        #prob_dict = Dict{SubString{String}, Dict{SubString{String}, Float64}}()
        prob_dict = Dict()
        for cls in states
            prob_dict[cls] = Dict()
        end
        for line in eachline(f)
            str_lin = split(strip(replace(line, r" +" => " ")))
            cur_c = str_lin[1]
            str_lin = parse.(Float64, str_lin[2:end]) 
            for (cls, prb) in zip(cols, str_lin)
                prob_dict[cur_c][cls] = prb
            end
        end
        return str_c, read, states, prob_dict
    end

    res = 1.
    for (s, c) in zip(str_c, read)
        res *= prob_dict[string(c)][string(s)]
    end

    open("output.txt", "w") do f
        println(f, string(res))
    end
end

end
