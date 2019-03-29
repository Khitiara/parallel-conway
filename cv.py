#!/usr/bin/python3
"""
Converts an alive count output file into a human-readable representation.
To plot this file in gnuplot, use the following commands:
> set view map
> set xrange [0:1024]
> set yrange [0:1024]
> splot "<filename>.dat" matrix with image
"""
import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print("Usage: cv <filename> <rowlen>")
    exit()

file_name = sys.argv[1]
rowlen = int(sys.argv[2])

data = [[0 for i in range(rowlen)] for j in range(rowlen)]

with open(file_name, 'rb') as file:
    for i in range(rowlen):
        line = file.read(rowlen * 4)
        for j in range(0, len(line), 4):
            data[i][j // 4] = int.from_bytes(line[j:j + 4], byteorder='little')

plt.imshow(data, cmap='hot', interpolation='nearest')
plt.savefig(file_name + '.png', frameon=False, dpi=640)