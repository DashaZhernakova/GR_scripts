#!/usr/bin/python
import sys

out = open(sys.argv[2], "w")
out.write("position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n")
with open(sys.argv[1]) as f:
    lines = [l.strip().split() for l in f.readlines()]
    for i in range(2, len(lines)):
        prev = lines[i-1]
        cur = lines[i]
        rate = (float(cur[2]) - float(prev[2])) / ((float(cur[3]) - float(prev[3]))/1000000)
        out.write(cur[3] + " " + str(rate) + " " + cur[2] + "\n")

out.close()