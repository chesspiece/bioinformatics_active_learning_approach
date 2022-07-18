import numpy as np

from five.dp_lib import dp_change

if __name__ == "__main__":
    with open("input_5_1.txt") as f:
        money = int(f.readline().strip())
        coins = np.array([int(x) for x in f.readline().strip().split()], dtype=np.int64)
    answ = dp_change(money, coins)

    with open("output.txt", "w") as f:
        f.write(" ".join(str(answ)))
