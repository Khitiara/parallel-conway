#!/usr/bin/python3
"""
Converts a Life output file into a human-readable representation.
"""
import sys

if len(sys.argv) != 3:
    print("Usage: cv <filename> <rowlen>")
    exit()

file_name = sys.argv[1]
rowlen = int(sys.argv[2])

with open(file_name, 'r') as file:
    with open(file_name + '.out', 'w') as out_file:
        for i in range(rowlen):
            line = file.read(rowlen)
            line = ''.join(map(lambda c: '.' if c == '\0' else '0', line))
            out_file.write(line)
            out_file.write('\n')
