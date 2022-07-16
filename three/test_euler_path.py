from copy import deepcopy

from three.grph_lib import Euler_grph, find_euler_path

if __name__ == "__main__":
    res = Euler_grph("input_3_8.txt")
    answ = find_euler_path(deepcopy(res.edges_list), res.i_deg, res.o_deg, res.vertex_list)
    with open("output.txt", "w") as f:
        f.write(" ".join(answ))
