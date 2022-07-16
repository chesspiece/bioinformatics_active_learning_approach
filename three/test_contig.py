import sys

from copy import deepcopy

from three.grph_lib import deBrujin_grph, max_nonbranch_paths, helper_dna_path

if __name__ == "__main__":
    sys.setrecursionlimit(30000)
    with open("input_3_12.txt", "r") as f:
        dna_list = f.readline().strip().split()

    res = deBrujin_grph(dna_list)
    # answ = find_euler_path(
    #    deepcopy(res.edges_list), res.i_deg, res.o_deg, res.vertex_list
    # )
    answ_tmp = max_nonbranch_paths(res.edges_list, res.vertex_list, res.i_deg, res.o_deg)
    answ = []
    for i in answ_tmp:
        answ.append(helper_dna_path(i))
    with open("output.txt", "w") as f:
        f.write(" ".join(answ))
