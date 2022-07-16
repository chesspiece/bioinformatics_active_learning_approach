from copy import deepcopy

from one.common_functions import neighbours

from three.grph_lib import deBrujin_grph, find_euler_path

if __name__ == "__main__":
    with open("input_3_10.txt", "r") as f:
        k = int(f.readline().strip())

    bin_list = neighbours("0"*k, k, choice_pat=["0", "1"])
    res = deBrujin_grph(bin_list)
    #answ = find_euler_path(
    #    deepcopy(res.edges_list), res.i_deg, res.o_deg, res.vertex_list
    #)
    answ = res.find_dna_cycle()
    with open("output.txt", "w") as f:
        f.write("".join(answ))
