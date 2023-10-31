from itertools import permutations
import subprocess
from random import randrange
i = 0
xs = " ".join(str(i) for i in range(5))
for perm in range(1000):
    xs = tuple(randrange(1, 1000) for _ in range(randrange(1, 100)))
    with open("t.txt", "w") as f:f.write(" ".join(str(i) for i in xs))
    result = subprocess.run(["./a.out", "t.txt"], capture_output=True, text=True)
    if " ".join(str(i) for i in sorted(xs)) != result.stdout.strip():
        # for i in range(len(xs)):
        #     if xs[i] != result.stdout[i]:
        #         print(i)
        #         break
        print("wtf", xs)
        # print(result.stdout)
        # print(xs)
        break
    i += 1
    if i%10 == 0:
        print(i)

else:
    print("Success!")
