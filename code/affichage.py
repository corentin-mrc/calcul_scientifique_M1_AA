import matplotlib.pyplot as plt

N = 25
I = (N - 1)**2
filename = "build/solution_exacte.txt"
with open(filename) as f:
    content = f.readlines()
    wk = []
    for line in content:
        if line != '' and line != '\n':
            wk.append(float(line))

