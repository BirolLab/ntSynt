#!/usr/bin/env python3
"""
extract_synteny_regions.py

For each genome assembly in a multi-genome ntSynt synteny block TSV, extract:
  1. Syntenic regions         -> <genome_basename>.syntenic.fa
  2. Non-syntenic regions     -> <genome_basename>.non_syntenic.fa

Dependencies: bedtools, samtools (for faidx), Python 3.6+

Usage:
    python extract_synteny_regions.py \
        --synteny  synteny_blocks.tsv \
        --fastas   genome1.fa genome2.fa genome3.fa \
        --outdir   output/

The --fastas list is used to resolve genome file names stored in the TSV
(column 2) to actual paths on disk.  Matching is done by basename, so
genome1.fa and /some/path/genome1.fa are treated as the same file.
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: str, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and optionally raise on failure."""
    print(f"  [cmd] {cmd}", flush=True)
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            text=True,
            stderr=subprocess.PIPE,
            check=check,
        )
    except subprocess.CalledProcessError as e:
        # When check=True subprocess raises CalledProcessError; present stderr
        print(f"ERROR: command failed (exit {e.returncode}):\n{e.stderr}", file=sys.stderr)
        sys.exit(1)
    return result


def fai_path(fasta: str) -> str:
    """Given a FASTA path, return the .fai index path."""
    return fasta + ".fai"


def ensure_fai(fasta: str) -> None:
    """Create a samtools .fai index if one doesn't already exist."""
    if not os.path.exists(fai_path(fasta)):
        print(f"  Indexing {fasta} with samtools faidx ...")
        run(f"samtools faidx {fasta}")


# ---------------------------------------------------------------------------
# Parse ntSynt TSV
# ---------------------------------------------------------------------------

def parse_synteny_tsv(tsv_path: str) -> dict:
    """
    Parse the ntSynt-style synteny blocks TSV.

    Columns (1-based):
      1  block_id
      2  genome_filename
      3  chrom
      4  start  (0-based; ntSynt uses BED-like coordinates)
      5  end
      6  strand
      7  n_minimizers
      8  discontinuity_reason

    Returns:
        { genome_filename: [(chrom, start, end), ...] }
    """
    blocks = defaultdict(list)

    with open(tsv_path, 'r', encoding='utf-8') as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                print(f"  WARNING: skipping malformed line {lineno}: {line}",
                      file=sys.stderr)
                continue
            genome = parts[1]
            chrom  = parts[2]
            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                print(f"  WARNING: non-integer coordinates on line {lineno}, skipping.",
                      file=sys.stderr)
                continue
            blocks[genome].append((chrom, start, end))

    return dict(blocks)


# ---------------------------------------------------------------------------
# Write BED
# ---------------------------------------------------------------------------

def write_bed(intervals, bed_path: str, genome_file: str) -> None:
    """Write intervals to a BED file. Intervals are assumed non-overlapping."""
    with open(f"{bed_path}.tmp", "w", encoding='utf-8') as fh:
        for chrom, start, end in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\n")
    run(f"bedtools sort -i {bed_path}.tmp -g {genome_file} > {bed_path}")
    run(f"rm {bed_path}.tmp")


# ---------------------------------------------------------------------------
# Extract FASTA sequences
# ---------------------------------------------------------------------------

def getfasta(fasta: str, bed: str, out_fa: str) -> None:
    """Run bedtools getfasta to pull sequences defined by a BED file."""
    run(f"bedtools getfasta -fi {fasta} -bed {bed} -fo {out_fa}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Main - extract syntenic and non-syntenic regions for each genome in the TSV."""
    parser = argparse.ArgumentParser(
        description="Extract syntenic and non-syntenic regions per genome assembly."
    )
    parser.add_argument(
        "--synteny", required=True,
        help="ntSynt-style synteny blocks TSV file."
    )
    parser.add_argument(
        "--fastas", required=True, nargs="+",
        help="One or more genome FASTA files. Basenames must match column 2 of the TSV."
    )
    parser.add_argument(
        "--outdir", default=".",
        help="Directory for output files (default: current directory)."
    )
    parser.add_argument(
        "--keep-bed", action="store_true",
        help="Keep intermediate BED files (useful for debugging)."
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Map basename -> full path for every supplied FASTA
    fasta_map = {}
    for fa in args.fastas:
        fasta_map[os.path.basename(fa)] = fa

    # Parse synteny TSV
    print("Parsing synteny TSV ...")
    synteny = parse_synteny_tsv(args.synteny)

    missing = set(synteny.keys()) - set(fasta_map.keys())
    if missing:
        print(
            "WARNING: the following genomes appear in the TSV but no matching FASTA "
            f"was supplied:\n  {', '.join(sorted(missing))}",
            file=sys.stderr
        )

    # Process each genome
    for genome_name, intervals in synteny.items():
        if genome_name not in fasta_map:
            print(f"  Skipping {genome_name} (no FASTA provided).")
            continue

        fasta = fasta_map[genome_name]
        base  = os.path.splitext(genome_name)[0]   # strip .fa / .fasta
        print(f"\n{'='*60}")
        print(f"Processing {genome_name}  ({len(intervals)} syntenic intervals)")
        print(f"{'='*60}")

        # 1. Index FASTA (needed by both bedtools complement and getfasta)
        ensure_fai(fasta)

        # 2. Write syntenic BED directly from TSV intervals (non-overlapping)
        syntenic_bed = os.path.join(args.outdir, f"{base}.syntenic.bed")
        print(f"  Writing syntenic BED -> {syntenic_bed}")
        write_bed(intervals, syntenic_bed, fai_path(fasta))

        # 3. Non-syntenic BED via bedtools complement
        non_syntenic_bed = os.path.join(args.outdir, f"{base}.non_syntenic.bed")
        print(f"  Computing non-syntenic BED (complement) -> {non_syntenic_bed}")
        run(f"bedtools complement -i {syntenic_bed} -g {fai_path(fasta)} > {non_syntenic_bed}")

        # 4. Extract FASTAs
        syntenic_fa     = os.path.join(args.outdir, f"{base}.syntenic.fa")
        non_syntenic_fa = os.path.join(args.outdir, f"{base}.non_syntenic.fa")

        print(f"  Extracting syntenic sequences     -> {syntenic_fa}")
        getfasta(fasta, syntenic_bed, syntenic_fa)

        print(f"  Extracting non-syntenic sequences -> {non_syntenic_fa}")
        getfasta(fasta, non_syntenic_bed, non_syntenic_fa)

        # 5. Optionally remove intermediate BED files
        if not args.keep_bed:
            for f in (syntenic_bed, non_syntenic_bed):
                os.remove(f)

    print("\nDone.")


if __name__ == "__main__":
    main()
