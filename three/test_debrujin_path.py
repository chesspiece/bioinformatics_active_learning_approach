from copy import deepcopy

from three.grph_lib import deBrujin_grph, find_euler_path

if __name__ == "__main__":
    with open("input_3_9.txt", "r") as f:
        k = int(f.readline().strip())
        dna_list = f.readline().strip().split()

    res = deBrujin_grph(dna_list)
    #answ = find_euler_path(
    #    deepcopy(res.edges_list), res.i_deg, res.o_deg, res.vertex_list
    #)
    answ = res.find_dna()
    with open("output.txt", "w") as f:
        f.write("".join(answ))
