from matplotlib import pyplot
import json
from math import log2
data = {}

with open("points for 4.json", "r") as f:
    data = json.load(f)

pyplot.xlabel("Number of planets (N)", fontsize=20)
pyplot.ylabel("Time taken in ms (T)",fontsize=20)
pyplot.title("N vs T graph for 90000 queries", fontsize=30)
for (x, y) in data["90000"]:
    pyplot.plot( x, y, "o", color='0.05')
# for (x, y) in data["20000"]:
#     pyplot.plot( x, y, "o", color='0.1')
# for (x, y) in data["30000"]:
#     pyplot.plot( x, y, "o", color='0.15')
# for (x, y) in data["40000"]:
#     pyplot.plot( x, y, "o", color='0.20')
# for (x, y) in data["50000"]:
#     pyplot.plot( x, y, "o", color='0.25')
# for (x, y) in data["60000"]:
#     pyplot.plot( x, y, "o", color='0.30')
#     # pyplot.plot( x, log2(x*1/5)*6, "ro")
# for (x, y) in data["70000"]:
#     pyplot.plot( x, y, "o", color='0.35')
# for (x, y) in data["80000"]:
#     pyplot.plot( x, y, "o", color='0.40')
# for (x, y) in data["90000"]:
#     pyplot.plot( x, y, "o", color='0.45')
    # pyplot.plot( x, log2(x*1/5)*6, "ro")

# pyplot.xscale('log')
pyplot.show()