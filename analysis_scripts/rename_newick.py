#!/usr/bin/env python3
'''
Rename a newick file based on the provided TSV
first column: old name
second column: new name
'''
from functools import reduce
import sys

def read_rename_tsv(rename_tsv):
    "Read the rename TSV into a dictionary"
    rename_dict = {}
    with open(rename_tsv, 'r', encoding="utf-8") as fin:
        for line in fin:
            old_name, new_name = line.strip().split('\t')
            rename_dict[old_name] = new_name

    return rename_dict

def rename_newick(newick, rename_dict):
    "Rename the newick file based on the dictionary"
    with open(newick, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip()
            line = reduce(lambda s, kv: s.replace(*kv), rename_dict.items(), line)
            print(line)

def main():
    "Rename the given newick file"
    if len(sys.argv[1:]) != 2:
        print(f"Usage: {sys.argv[0]} <newick file (- if read from stdin)> <rename TSV>")
        sys.exit()

    newick = "/dev/stdin" if sys.argv[1] == '-' else sys.argv[1]
    rename_tsv = sys.argv[2]

    rename_dict = read_rename_tsv(rename_tsv)

    rename_newick(newick, rename_dict)

if __name__ == '__main__':
    main()
