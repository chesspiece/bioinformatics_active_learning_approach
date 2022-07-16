import sys

from three.grph_lib import deBrujin_grph2

if __name__ == "__main__":
    sys.setrecursionlimit(30000)
    with open("input_3_11.txt", "r") as f:
        k, d = [int(x) for x in f.readline().strip().split()]
        dna_list = f.readline().strip().split()

    res = deBrujin_grph2(dna_list)
    # answ = find_euler_path(
    #    deepcopy(res.edges_list), res.i_deg, res.o_deg, res.vertex_list
    # )
    answ = res.find_dna(d)
    with open("output.txt", "w") as f:
        f.write("".join(answ))
