from copy import deepcopy

from three.grph_lib import Euler_grph, find_euler_cycle

if __name__ == "__main__":
    res = Euler_grph("input_3_6.txt")
    answ = find_euler_cycle(deepcopy(res.edges_list), res.vertex_list[0])
    with open("output.txt", "w") as f:
        f.write(" ".join(answ))
