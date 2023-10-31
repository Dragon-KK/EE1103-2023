import subprocess
import sys


with open("outputs/" + sys.argv[2] + "_output.txt", "r") as f:
    subprocess.run(("a" + sys.argv[1] + ".exe", "inputs/" + sys.argv[2] + ".txt"))

    with open("ee23b135_quiz2_q" + sys.argv[1] + "_output.txt") as f2:
        xs = f2.read()
    xs = xs.split("\n")
    if xs[0].strip() == f.readline().strip():
        print("Success!", xs[1])
    else:
        print("Fail", xs[1])
        # print(xs[0])