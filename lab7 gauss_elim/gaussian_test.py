import numpy as np

import subprocess

import random
random.seed(19)
for _ in range(101, 400):
    n = _//100
    A = np.array([
        [random.randrange(0, 2000)/10 for _ in range(n)] for __ in range(n)
    ])
    B = np.array(
        [random.randrange(0, 10000)/10 for _ in range(n)]
    )
    with open("t.dat", "w") as f:
        for i in range(n):
            f.write(" ".join(str(k) for k in (A[i])) + " " + str(B[i]) + "\n")
    X2 = np.linalg.solve(A,B)
    
    xs = subprocess.run(("a" + ".exe", "t.dat", str(n)), text=True, capture_output=True)
    arr = np.array([float(i) for i in xs.stdout.strip().split(" ")])
    
    if not all(np.isclose(X2, arr, rtol=0.0001, atol=0.00001, equal_nan=False)):
        print("VROROOO!", _)
        break
    print(_)
else:
    print("Success!")
