#!/usr/bin/env python3
"""
run_pipeline.py  –  driver for the multi-genome synteny Snakemake pipeline.

Builds the pipeline config from command-line arguments and launches Snakemake
programmatically, so no config.yaml needs to be maintained by hand.

Usage examples
--------------
# Dry-run (print what would be done)
python run_pipeline.py \
    --accessions assemblies.tsv \
    --group ichneumonidae \
    --date 250528 \
    --block-stats-script /path/to/ntSynt/analysis_scripts/denovo_synteny_block_stats.py \
    --mx-stats-script    /path/to/synteny/scripts/analyze_mx.py \
    --dry-run

# Real run, 16 cores
python run_pipeline.py \
    --accessions assemblies.tsv \
    --group ichneumonidae \
    --date 250528 \
    --block-stats-script /path/to/denovo_synteny_block_stats.py \
    --mx-stats-script    /path/to/analyze_mx.py \
    --cores 16

# With an optional tree file for ntSynt-viz
python run_pipeline.py ... --tree species.nwk

# Force-rerun specific rules
python run_pipeline.py ... --forcerun run_ntsynt plot_divergences
"""

import argparse
import datetime
import subprocess
import sys
import json
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

SNAKEFILE = Path(__file__).parent / "Snakefile"
SCRIPTS_DIR = Path(__file__).parent


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    p = argparse.ArgumentParser(
        description="Run the multi-genome synteny analysis pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ------------------------------------------------------------------
    # Main inputs
    # ------------------------------------------------------------------
    req = p.add_argument_group("main inputs")
    req.add_argument(
        "--accessions", required=False, metavar="TSV",
        help=(
            "TSV listing NCBI genome accessions to use - one per line"
        ),
    )
    req.add_argument(
        "--genomes", required=False, metavar="TSV",
        help=(
            "TSV listing paths to genome assemblies to analyze - one per line"
        )
    )
    req.add_argument(
        "--prefix", required=True, metavar="NAME",
        help="Prefix name for genome assemblies (e.g. 'ichneumonidae').",
    )
    req.add_argument(
        "--name-conversions", required=True, metavar="CONVERSIONS",
        help=(
            "TSV file listing name conversions for display purposes. "
            "Expected columns: accession/assembly base name; new name"
        )
    )

    # ------------------------------------------------------------------
    # Optional ntSynt / analysis parameters
    # ------------------------------------------------------------------
    opt = p.add_argument_group("analysis parameters")
    opt.add_argument("--fpr",  type=float, default=0.025,
                     help="Bloom filter false-positive rate for ntSynt.")
    opt.add_argument("--ntsynt-k", type=int, default=24,
                     help="k-mer size for ntSynt.")
    opt.add_argument("--ntsynt-w", type=int, default=1000,
                     help="Minimizer window size for ntSynt.")
    opt.add_argument("--make-tree", action="store_true",
                     help=(
                         "Automatically generate phylogenetic tree for ribbon plot. "
                         "If --accessions specified, will look for mitochondrial sequences in these downloaded files. "
                         "If --mt-genomes supplied, will use those mitochondrial genomes. "
                         "Otherwise, will use nuclear genomes with Mash + Quicktree. "
                     ))
    opt.add_argument("--mt-genomes", required=False, metavar="FASTA",
                     help=(
                         "FASTA file containing mitochondrial genomes for all input assemblies. "
                         "If specified, must also supply a TSV file for converting header names to new names."
                         ))
    opt.add_argument("--mt-name-conversions", required=False, metavar="TSV",
                     help=(
                         "TSV listing name conversions between mitochondrial genome accessions and assembly names. "
                         "Expected columns: mt genome accession; new name (matching --name-conversions)"
                     ))
    opt.add_argument("--tree", default="", metavar="NEWICK",
                     help="Optional Newick tree file for ntSynt-viz. Omit to skip.")
    opt.add_argument("--ntsynt-viz_ribbon-adjust", type=float, default=0.2,
                     help="Adjustment factor for ntSynt-viz ribbons. Increase if ribbon plot labels are cut off.")
    # ------------------------------------------------------------------
    # Snakemake execution options
    # ------------------------------------------------------------------
    smk = p.add_argument_group("snakemake options")
    smk.add_argument("--cores", type=int, default=12,
                     help="Number of CPU cores to use.")
    smk.add_argument("--dry-run", "-n", action="store_true",
                     help="Perform a dry run (print rules, do not execute).")
    smk.add_argument("--forcerun", nargs="*", metavar="RULE",
                     help="Force re-execution of specific rules (space-separated).")
    smk.add_argument("--until", nargs="*", metavar="RULE",
                     help="Run the pipeline only up to and including these rules.")
    smk.add_argument("--snakemake-args", nargs=argparse.REMAINDER,
                     default=[], metavar="...",
                     help=(
                         "Any additional arguments passed verbatim to Snakemake "
                         "(place after all other flags, e.g. -- --rerun-incomplete)."
                     ))

    return p.parse_args(), p


def build_config(args: argparse.Namespace) -> dict:
    """Translate parsed CLI args into the config dict the Snakefile expects."""
    return {
        "accessions":    str(Path(args.accessions).resolve()) if args.accessions else "",
        "genomes": str(Path(args.accessions).resolve()) if args.genomes else "",
        "prefix":  args.prefix,
        "name_conversions": args.name_conversions,
        "make_tree": args.make_tree,
        "mt_genomes": args.mt_genomes,
        "mt_name_conversions": args.mt_name_conversions,
        "date":              f"{args.prefix}_assemblies",
        "fpr":               args.fpr,
        "ntsynt_k":          args.ntsynt_k,
        "ntsynt_w":          args.ntsynt_w,
        "treefile":          args.tree,
        "ntsynt_viz_ribbon_adjust": args.ntsynt_viz_ribbon_adjust,
        "scripts_dir":     str(SCRIPTS_DIR.resolve()),
    }


def validate_paths(args: argparse.Namespace) -> None:
    """Abort early if required input files are missing."""
    errors = []
    if args.genomes and not Path(args.genomes).exists():
        errors.append(f" --genomes: file not found: {args.genomes}")
    elif args.genomes:
        with open(args.genomes, 'r', encoding="utf-8") as fin:
            for genome in fin:
                if not Path(genome).exists():
                    errors.append(f"Genome file listed in --genome not found: {genome}")
    if args.accessions and not Path(args.accessions).exists():
        errors.append(f" --accessions: file not found: {args.accessions}")
    if args.tree and not Path(args.tree).exists():
        errors.append(f"  --tree: file not found: {args.tree}")
    if args.make_tree and args.mt_genomes and not Path(args.mt_genomes).exists():
        errors.append(f"  --mt-genomes: file not found: {args.mt_genomes}")
    if args.make_tree and args.mt_genomes and not Path(args.mt_name_conversions).exists():
        errors.append(f"  --mt-name-conversions: file not found: {args.mt_name_conversions}")
    if errors:
        print("ERROR: the following required files were not found:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_snakemake_cmd(args: argparse.Namespace, config: dict) -> list[str]:
    """Construct the Snakemake command-line invocation as a list of tokens."""
    cmd = [
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--cores",     str(args.cores),
        "--printshellcmds",
        "--nolock",
    ]

    # One --config flag followed by all key=value pairs as separate tokens
    config_pairs = []
    for k, v in config.items():
        if isinstance(v, dict):
            config_pairs.append(f"{k}={json.dumps(v)}")
        else:
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

def validate_options(args, parser):
    """Validate that input arguments are compatible"""
    if not args.accessions and not args.genomes:
        raise parser.error("Please specify either --accessions or --genomes")
    if args.accessions and args.genomes:
        raise parser.error("Please specify one of --accessions or --genomes")
    if args.mt_genomes and not args.mt_name_conversions:
        raise parser.error("If --mt-genomes is supplied, please also supply --mt-name-conversions")
    if not args.mt_genomes and args.mt_name_conversions:
        print("WARNING: --mt-name-conversions only used when --mt-genomes specified")
    if not args.make_tree and args.mt_genomes:
        print("WARNING: --mt-genomes is only used when --make-tree is specified.")
    if args.tree and args.make_tree:
        print("WARNING: --tree specified, so will override --make-tree.")
        args.make_tree = False

def main() -> None:
    """Main entry point: parse args, validate, build config, and launch Snakemake."""
    args, parser = parse_args()
    validate_options(args, parser)
    validate_paths(args)
    config = build_config(args)

    print("=" * 60)
    print("Synteny pipeline  –  effective configuration")
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

    result = subprocess.run(cmd, check=False)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
