#!/usr/bin/env python3
'''
Analyze the minimizer density, and minimizer overlap for an ntSynt run
'''
import argparse
from datetime import datetime
import re
import sys
import subprocess
import os
from typing import Iterator, Set
import tempfile
import btllib


def count_minimizers_graph(dot_file):
    "Count minimizers in DOT graph"
    node_pattern = re.compile(r'^"\d+"\s+\[label=')
    count_mx = 0
    with open(dot_file, 'r', encoding="utf-8") as f:
        for line in f:
            if node_pattern.match(line):
                count_mx += 1
    return count_mx

def reverse_complement_kmer(kmer):
    """Return canonical representation of a k-mer."""
    trans = str.maketrans("ACGT", "TGCA")
    return kmer.translate(trans)[::-1]


def find_shared_kmers(fastas, k):
    """Find canonical k-mers shared between input FASTA files using KMC.

    Args:
        fastas (str): Path to file with fasta file paths - 1 per line.
        k (int): k-mer size.

    Returns:
        set[str]: Set of common k-mers between the input fastas
    """
    if k <= 0:
        raise ValueError("k must be a positive integer")

    shared_kmers = set()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_dir, "kmc_find_common_kmers.sh")

    # Use a temporary directory to avoid clutter
    with tempfile.TemporaryDirectory(dir=os.getcwd(), prefix="kmc_") as tmp_dir:
        out_prefix = os.path.join(tmp_dir, "shared_kmers_kmc")

        cmd = [
            script_path,
            fastas,
            str(k),
            out_prefix,
            tmp_dir,
        ]

        # Run the KMC pipeline
        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"KMC pipeline failed:\nSTDOUT:\n{exc.stdout}\nSTDERR:\n{exc.stderr}"
            ) from exc

        # Read output k-mers
        output_file = f"{out_prefix}.kmers.txt"

        if not os.path.exists(output_file):
            raise RuntimeError("Expected output file not found after KMC run")

        with open(output_file, "r", encoding="utf-8") as handle:
            for line in handle:
                kmer = line.strip().split("\t")[0]
                if kmer:
                    shared_kmers.add(kmer)

    return shared_kmers


def iter_minimizers(path: str) -> Iterator[str]:
    """Yield minimizer hashes from a sketch file."""
    with open(path, "r", encoding="utf-8") as mx_in:
        for line in mx_in:
            parts = line.strip("\n").split("\t")
            if len(parts) < 2:
                continue

            for entry in parts[1].split():
                yield entry.split(":")[0]


def find_shared_minimizers(file) -> Set[str]:
    """Find minimizers shared across multiple sketch files."""
    with open(file, 'r', encoding="utf-8") as f:
        file_list = [line.strip() for line in f if line.strip()]
    if not file_list:
        return set()

    # Initialize with first file
    common = set(iter_minimizers(file_list[0]))

    # Intersect with remaining files
    for path in file_list[1:]:
        current = {
            m for m in iter_minimizers(path)
            if m in common
        }

        common = current

        if not common:
            break

    return common

def count_minimizers_sketch(mx_file, common_bf, common_hash):
    "Compute the number of mx in the sketches that are in the BF and common hash"
    num_mx, num_in_bf, num_in_hash = 0, 0, 0
    with open(mx_file, 'r', encoding="utf-8") as mx_in:
        for line in mx_in:
            parts = line.strip().split("\t")
            entries = parts[1].split()

            for entry in entries:
                seq = entry.split(":")[-1]
                seq = seq.upper()
                num_mx += 1

                if common_bf.contains(seq):
                    num_in_bf += 1
                if seq in common_hash or reverse_complement_kmer(seq) in common_hash:
                    num_in_hash += 1
    return num_mx, num_in_bf, num_in_hash

def load_shared_kmers(filepath):
    """Load shared k-mers from a file into a set."""
    shared_kmers = set()
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            shared_kmers.add(line.strip().split("\t")[0])
    return shared_kmers

def main():
    "Analyze minimizer density and overlap"
    parser = argparse.ArgumentParser(description="Analyze minimizer density and overlap")
    parser.add_argument("--mx", help="Path to file with minimizer files - one per line", required=True, type=str)
    parser.add_argument("-k", help="k-mer size", required=True, type=int, default=24)
    parser.add_argument("--bf", help="Path to common BF", required=True, type=str)
    parser.add_argument("--fastas", help="Path to file with FASTA files - one per line", required=True, type=str)
    parser.add_argument("--dot", help="Path to graph file", required=True, type=str)
    parser.add_argument("--shared_kmers", help="Path to output file for shared k-mers", required=False, type=str)
    args = parser.parse_args()

    print(f"{datetime.now()} Analyzing dot graph..", file=sys.stderr)
    num_mx_dot = count_minimizers_graph(args.dot)

    print(f"{datetime.now()} Loading common Bloom filter..", file=sys.stderr)
    common_bf = btllib.KmerBloomFilter(args.bf)
    print(f"{datetime.now()} Finding shared k-mers..", file=sys.stderr)
    common_hash = load_shared_kmers(args.shared_kmers) if args.shared_kmers else find_shared_kmers(args.fastas, args.k)

    print(f"{datetime.now()} Finding shared minimizers..", file=sys.stderr)
    common_minimizers = find_shared_minimizers(args.mx)

    print(f"{datetime.now()} Counting minimizers in sketches..", file=sys.stderr)
    results = []
    with open(args.mx, 'r', encoding="utf-8") as f:
        for line in f:
            results.append(count_minimizers_sketch(line.strip(), common_bf, common_hash))

    min_entry = min(results, key=lambda x: x[0])
    max_entry = max(results, key=lambda x: x[0])

    print("num_mx_in_dot", "common_mx_sketch", "max_num_mx", "max_num_mx_bf", "max_num_mx_hash",
          "min_num_mx", "min_num_mx_bf", "min_num_mx_hash", sep="\t")
    print(num_mx_dot, len(common_minimizers),
          max_entry[0], max_entry[1], max_entry[2],
          min_entry[0], min_entry[1], min_entry[2], sep="\t")



if __name__ == "__main__":
    main()
