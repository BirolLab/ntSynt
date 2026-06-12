#!/usr/bin/env python3
"""
run_phylogeny.py  —  driver for the phylogenetic tree Snakemake pipeline.

Builds config from command-line arguments, auto-fills paths shared with the
master pipeline, and launches phylogeny.smk.

Usage examples
--------------
# Auto-detect MT source - download more MT if needed, or fall-back on nuclear otherwise
python run_phylogeny.py \
    --group lucinidae \
    --date  26-06-02 

# Force nuclear tree
python run_phylogeny.py \
    --group lucinidae --date 26-06-02 \
    --mt-source nuclear 

# User-supplied MT FASTA
python run_phylogeny.py \
    --group lucinidae --date 26-06-02 \
    --mt-source user-fasta \
    --mt-fasta  /path/to/mt_assemblies.fa 
"""

import argparse
from email import errors
import json
import os
import subprocess
import sys
from pathlib import Path

SNAKEFILE   = Path(__file__).parent / "phylogeny.smk"
SCRIPTS_DIR = Path(__file__).parent


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run the phylogenetic tree pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ------------------------------------------------------------------
    # Shared with master pipeline
    # ------------------------------------------------------------------
    req = p.add_argument_group("required inputs")
    req.add_argument("--group", required=True, metavar="NAME",
                     help="Taxonomic group name (must match master pipeline run).")
    req.add_argument("--date", required=True, metavar="YY-MM-DD",
                     help="Date string used in the master pipeline download directory.")

    # ------------------------------------------------------------------
    # MT source control
    # ------------------------------------------------------------------
    mt = p.add_argument_group("MT source")
    mt.add_argument(
        "--mt-source", default="auto",
        choices=["auto", "user-fasta", "nuclear"],
        help=(
            "MT sequence source. "
            "'auto': detect from sequence report, then try to download based on species, fall back to nuclear. "
            "'user-fasta': use --mt-fasta. "
            "'nuclear': skip MT, use mashtree on nuclear genomes."
        ),
    )
    mt.add_argument("--mt-fasta", default="", metavar="PATH",
                    help="User-supplied MT FASTA (required when --mt-source user-fasta).")

    # ------------------------------------------------------------------
    # Optional overrides
    # ------------------------------------------------------------------
    opt = p.add_argument_group("optional parameters")
    opt.add_argument("--prefix", default=None,
                     help="Output file prefix. Defaults to --group.")
    opt.add_argument("--threads", type=int, default=12,
                     help="Number of threads for mafft, iqtree, mashtree.")

    # ------------------------------------------------------------------
    # Snakemake options
    # ------------------------------------------------------------------
    smk = p.add_argument_group("snakemake options")
    smk.add_argument("--cores", type=int, default=12)
    smk.add_argument("--dry-run", "-n", action="store_true")
    smk.add_argument("--forcerun", nargs="*", metavar="RULE")
    smk.add_argument("--until", nargs="*", metavar="RULE")
    smk.add_argument("--snakemake-args", nargs=argparse.REMAINDER, default=[])

    return p.parse_args()


def build_config(args: argparse.Namespace) -> dict:
    fam_low      = args.group.lower()
    assembly_dir = f"{args.date}_assemblies"

    # Paths produced by the master pipeline
    seq_report   = f"{assembly_dir}/{fam_low}_sequence-reports.tsv"
    name_conv    = f"{assembly_dir}/{fam_low}_name_conversion.tsv"
    fasta_list   = f"{assembly_dir}/{fam_low}_fasta_list.txt"

    return {
        "taxonomic_group":      fam_low,
        "date":                 args.date,
        "seq_report":           str(Path(seq_report).resolve()),
        "name_conversion":      str(Path(name_conv).resolve()),
        "fasta_list":           str(Path(fasta_list).resolve()),
        "mt_source":            args.mt_source,
        "mt_fasta":             args.mt_fasta,
        "prefix":               args.prefix or fam_low,
        "threads":              args.threads,
        "scripts_dir":     str(SCRIPTS_DIR.resolve()),
    }


def validate_paths(args: argparse.Namespace) -> None:
    error = None
    if args.mt_source == "user-fasta":
        if not args.mt_fasta:
            error = "  --mt-fasta is required when --mt-source=user-fasta"
        elif not Path(args.mt_fasta).exists():
            error = f"  --mt-fasta: not found: {args.mt_fasta}"

    if error:
        print("ERROR:", file=sys.stderr)
        print(error, file=sys.stderr)
        sys.exit(1)


def build_snakemake_cmd(args: argparse.Namespace, config: dict) -> list[str]:
    cmd = [
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--cores",     str(args.cores),
        "--printshellcmds",
    ]

    config_pairs = []
    for k, v in config.items():
        config_pairs.append(f"{k}={v}")

    if config_pairs:
        cmd += ["--config"] + config_pairs

    if args.dry_run:
        cmd.append("--dry-run")

    for rule in (args.forcerun or []):
        cmd += ["--forcerun", rule]

    for rule in (args.until or []):
        cmd += ["--until", rule]

    cmd += args.snakemake_args
    return cmd


def main() -> None:
    args = parse_args()
    validate_paths(args)
    config = build_config(args)

    print("=" * 60)
    print("Phylogeny pipeline  —  effective configuration")
    print("=" * 60)
    for k, v in config.items():
        if isinstance(v, dict):
            for sk, sv in v.items():
                print(f"  scripts.{sk:<28} {sv}")
        else:
            print(f"  {k:<32} {v}")
    print("=" * 60)
    if args.dry_run:
        print("DRY RUN — no files will be created.\n")

    cmd = build_snakemake_cmd(args, config)
    print("Running:", " ".join(cmd), flush=True)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()