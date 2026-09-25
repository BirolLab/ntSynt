#!/usr/bin/env python3
'''
Given a set of synteny blocks, compute basic stats including the number of blocks,
coverage and average/median length
'''
import argparse
import re
import os
from collections import namedtuple
import numpy as np


Block = namedtuple("Block", ["length", "block_id"])

def get_lengths(list_block: list, block_tally=None, asm_threshold=0) -> list:
    "From a list of Block, get the lengths that satisfy the asm_threshold"
    if block_tally is None:
        return [block.length for block in list_block]
    return [block.length for block in list_block if block_tally[block.block_id] >= asm_threshold]

def read_blocks(tsv_file: str) -> tuple:
    "Given the synteny block TSV file, return a dict of average block lengths"
    block_lengths = {} # asm -> [block_length]
    block_id_tallies = {} # block_id -> set(assemblies)

    with open(tsv_file, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            block_id, asm, start, end = line[0], line[1], line[3], line[4]
            start, end = int(start), int(end)
            if asm not in block_lengths:
                block_lengths[asm] = []

            block_lengths[asm].append(Block(end - start, block_id))
            if block_id not in block_id_tallies:
                block_id_tallies[block_id] = set()
            block_id_tallies[block_id].add(asm)

    # Convert block ID tallies to block_id -> count distinct assemblies
    block_id_tallies = {block_id: len(set_asm) for block_id, set_asm in block_id_tallies.items()}

    return block_lengths, block_id_tallies

def calculate_ng50(blocks: list, genome_size: float) -> float:
    "Calculate the NG50 block length"
    target_size = genome_size *  0.5
    curr_sum = 0
    for block in sorted(blocks, reverse=True):
        curr_sum += block
        if curr_sum >= target_size:
            return block
    return 0

def get_genome_size(fai_filename: str) -> int:
    "Get the genome size from the input FAI file"
    curr_sum = 0
    with open(fai_filename, 'r', encoding="utf-8") as fin:
        for line in fin:
            length_ctg = int(line.strip().split("\t")[1])
            curr_sum += length_ctg
    return curr_sum


def read_fais(list_fais: list) -> dict:
    "Read the genome lengths per genome into a dictionary"
    fai_dict = {}
    re_fai_filename = re.compile(r'^(\S+).fai')
    for fai in list_fais:
        if re_match := re.search(re_fai_filename, fai):
            filename = os.path.basename(re_match.group(1))
            fai_dict[filename] = get_genome_size(fai)
    return fai_dict


def main() -> None:
    "Run the de novo stats assessment"
    parser = argparse.ArgumentParser(description="Compute de novo stats on synteny blocks")
    parser.add_argument("--tsv", help="ntSynt synteny block file", required=True)
    parser.add_argument("--fai", help="FAI files for the compared genomes", nargs="+",
                        required=True)
    args = parser.parse_args()

    genome_sizes = read_fais(args.fai)
    block_lengths, block_id_tallies = read_blocks(args.tsv)

    num_genomes = len(args.fai)
    num_blocks = sum((len(get_lengths(list_block)) for _, list_block in block_lengths.items()))/num_genomes
    num_blocks_all_asm = sum((len(get_lengths(list_block, block_id_tallies, asm_threshold=num_genomes)) \
        for _, list_block in block_lengths.items()))/num_genomes
    total_length = sum((sum(get_lengths(list_block)) for _, list_block in block_lengths.items()))/num_genomes
    avg_coverage = sum(
        (sum(get_lengths(list_block))/genome_sizes[asm]*100 for asm, list_block in block_lengths.items()))/num_genomes
    sorted_genome_sizes = sorted([(genome_size, asm) for asm, genome_size in genome_sizes.items()], key=lambda x:x[0])
    min_genome_size, min_genome_asm = sorted_genome_sizes[0]
    coverage_min_genome_size = (sum(get_lengths(block_lengths[min_genome_asm]))/min_genome_size) * 100
    avg_coverage_all = sum(
        (sum(get_lengths(list_block, block_id_tallies, asm_threshold=num_genomes))/genome_sizes[asm]*100 \
            for asm, list_block in block_lengths.items()))/num_genomes
    average_block_length = sum(
        (np.mean(get_lengths(list_block)) for _, list_block in block_lengths.items()))/num_genomes
    median_block_length = sum(
        (np.median(get_lengths(list_block)) for _, list_block in block_lengths.items()))/num_genomes
    average_ng50 = sum(
        (calculate_ng50(get_lengths(list_block), genome_sizes[asm]) \
            for asm, list_block in block_lengths.items()))/num_genomes
    average_n50 = sum(
        (calculate_ng50(get_lengths(list_block), sum(get_lengths(list_block))) \
            for _, list_block in block_lengths.items()))/num_genomes

    print("Number_genomes", "Number_blocks","Number_blocks_all_asm",  "Average_coverage", "Average_coverage_all_asm",
          "Coverage_min_genome_size", "Average_length", "Median_length", "Total_length",
          "NG50_length", "N50_length", sep="\t")
    print(f"{num_genomes}\t{int(num_blocks)}\t{int(num_blocks_all_asm)}\t{avg_coverage}\t{avg_coverage_all}\t"\
        f"{coverage_min_genome_size}\t{average_block_length}\t{median_block_length}\t{total_length}\t"\
        f"{int(average_ng50)}\t{int(average_n50)}")


if __name__ == "__main__":
    main()
