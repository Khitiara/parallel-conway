#!/usr/bin/python3
"""
Converts an alive count output file into a human-readable representation.
"""
import sys

if len(sys.argv) != 3:
    print("Usage: cv <filename> <rowlen>")
    exit()

file_name = sys.argv[1]
rowlen = int(sys.argv[2])

with open(file_name, 'rb') as file:
    with open(file_name + '.csv', 'w') as out_file:
        out_file.write('x,y,count\n');
        for i in range(rowlen):
            line = file.read(rowlen)
            for j in range(0, len(line), 4):
                count = int.from_bytes(line[j:j+4], byteorder='little')
                out_file.write('{0},{1},{2}\n'.format(j // 4, i, count));
