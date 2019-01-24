import random

num_datapoints = 500
num_vars = 100

from_file = open("./bnetflix.data", "r")
to_file = open("./netflix_{}.data".format(num_datapoints), "w+")

lines = []
for line in from_file:
    lines.append(line.replace(",", " "))

random.shuffle(lines)

to_file.write(str(num_vars) + "\n")
to_file.write("2 "*num_vars + "\n")
to_file.write(str(num_datapoints) + "\n")

for line in lines[:500]:
    to_file.write(line)