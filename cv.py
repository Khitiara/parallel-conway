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

if len(sys.argv) != 3:
    print("Usage: cv <filename> <rowlen>")
    exit()

file_name = sys.argv[1]
rowlen = int(sys.argv[2])

with open(file_name, 'rb') as file:
    with open(file_name + '.dat', 'w') as out_file:
        for i in range(rowlen):
            line = file.read(rowlen * 4)
            for j in range(0, len(line), 4):
                count = int.from_bytes(line[j:j+4], byteorder='little')
                out_file.write('{0} '.format(count))
            out_file.write('\n')
