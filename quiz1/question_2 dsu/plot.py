from random import randint
import subprocess
import sys

PROGRAM_NAME = "a" + sys.argv[1]

INPUT_TMP_FILE_NAME = PROGRAM_NAME + ".txt"
EXECUTABLE_NAME = PROGRAM_NAME + ".exe"
OUT_FILE_NAME = "ee23b135_quiz2_q" + PROGRAM_NAME[1] + "_output.txt"
DATA_FILE_NAME = PROGRAM_NAME + ".json"

def random_dat(n, m):
    return f"{n} {m}" + "\n" + "\n".join((f"{randint(0, n-1)} {randint(0, n-1)}") for _ in range(m))

datas = {}

m = 90000
datas[m] = []
n = 100
jump = 1
i = 0
while n <= 100000:
    i+=1
    with open(INPUT_TMP_FILE_NAME, "w") as f:
        f.write(random_dat(n, m))
    times_taken = []


    for _ in range(100):
        times_taken.append(int(subprocess.run((EXECUTABLE_NAME, INPUT_TMP_FILE_NAME), text=True, capture_output=True).stdout.split()[0]))

    times_taken.sort()

    datas[m].append((n, int(sum(times_taken[20:80])/60)))

    n += int(jump)
    if i//10:
        i = 0
        jump*=1.5
        
import json
with open(DATA_FILE_NAME, "w") as f:
    json.dump(datas, f)