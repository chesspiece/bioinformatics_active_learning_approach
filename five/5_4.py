import numpy as np

from five.dp_lib import Grph_path

if __name__ == "__main__":
    v_quant = 0
    e_quant = 0
    edges_list = []
    with open("input_5_4.txt") as f:
        source, sink = (int(x) for x in f.readline().strip().split())
        for line in f:
            inp, outp, weight = [int(x) for x in line.strip().split()]
            e_quant += 1
            v_quant = max(inp, outp, v_quant)
            edges_list.append([inp, outp, weight])
    v_quant += 1

    grph = Grph_path(v_quant, e_quant)
    grph.build_graph(edges_list)
    length, path = grph.longest_path(source, sink)

    with open("output.txt", "w") as f:
        f.write("".join(str(length) + "\n"))
        f.write(" ".join(map(str, path)))
