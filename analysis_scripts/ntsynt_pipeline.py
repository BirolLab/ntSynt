#!/usr/bin/env python3
"""
run_pipeline.py  –  driver for the multi-genome synteny Snakemake pipeline.

Builds the pipeline config from command-line arguments and launches Snakemake
programmatically.
"""

import argparse
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
        epilog=(
            "NOTE: All parameters available in ntSynt (https://github.com/BirolLab/ntsynt) "
            "and ntSynt-viz (https://github.com/BirolLab/ntSynt-viz) can also be supplied to this script"
        )
    )

    # ------------------------------------------------------------------
    # Main inputs
    # ------------------------------------------------------------------
    req = p.add_argument_group("main inputs")
    req.add_argument(
        "--accessions", required=False, metavar="TSV",
        help=(
            "TSV listing NCBI genome accessions to use - one per line.\n"
            "NOTE: Only chromosomes are retained after downloading accessions -"
            "if different behaviour is desired, use --genomes instead."
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
    opt.add_argument("--keep-incomplete", required=False, action="store_true",
                     help=(
                        "When --accessions is specified, by default the pipeline will skip analyzing any accessions "
                        "which appear incomplete (ie. there are unlocalized-scaffold entries with no assembled-molecule counterpart). "
                        "Use this option to still use these assemblies in the downstream analysis."
                    ))

    # ntSynt
    opt.add_argument("--hashes", type=int, default=None,
                    help="Number of hash functions for Bloom filter creation [ntSynt default: 3].")
    opt.add_argument("-b", "--block_size", help=argparse.SUPPRESS)
    opt.add_argument("--indel", help=argparse.SUPPRESS)
    opt.add_argument("--merge", help=argparse.SUPPRESS)
    opt.add_argument("--w_rounds", nargs="+", help=argparse.SUPPRESS, type=int)

    # ntSynt-viz
    opt.add_argument("--scale", help="Length of scale bar in bases for ntSynt-viz", type=float,
                     default=100e6)
    opt.add_argument("--seq_length", help="Minimum sequence length for ntSynt-viz", type=int)
    opt.add_argument("--ntsynt-viz_ribbon-adjust", default="auto",
                     help="Adjustment factor for ntSynt-viz ribbons. Increase if ribbon plot labels are cut off.")
    opt.add_argument("--target-genome", default="", metavar="NAME",
                    help="Target genome for ntSynt-viz (placed at top, ribbons coloured by its chromosomes).")
    opt.add_argument("--viz-length", type=int, default=None,
                    help="Minimum synteny block length for ntSynt-viz display (bp). Defaults to --block_size if set.")
    opt.add_argument("--format", choices=["png", "pdf", "svg"], default="png",
                    help="Output format for ntSynt-viz ribbon plot [png].")
    opt.add_argument("--width", type=float, default=None,
                    help="Width of ntSynt-viz ribbon plot in cm.")
    opt.add_argument("--dpi", type=int, default=None,
                    help="Resolution of ntSynt-viz ribbon plot (png only) [ntSynt-viz default: 300].")
    opt.add_argument("--centromeres", default="", metavar="TSV",
                    help="TSV file with centromere positions for ntSynt-viz (columns: bin_id, seq_id, start, end).")
    opt.add_argument("--order", default="", metavar="FILE",
                    help="File specifying genome order in ntSynt-viz ribbon plot.")
    opt.add_argument("--keep", nargs="+", default=None, metavar="GENOME:CHR",
                    help="Genome:chromosome pairs to show in ntSynt-viz (e.g. --keep genome1:chr1 genome2:chr3).")
    opt.add_argument("--no-arrow", action="store_true",
                    help="Do not draw strand-flip arrows in ntSynt-viz for normalization.")
    opt.add_argument("--optimize-ordering", action="store_true",
                    help="Optimize tree-guided genome sorting using inversions. "
                    "Only use with strictly bifurcating trees.")
    opt.add_argument("--haplotypes", type=str, default=None, metavar="haplotypes.tsv",
                    help="Optional TSV file listing haplotype information for each genome. If --genomes "
                    "is used, each row should be the two NEW assembly names (based on name conversion), separated by a tab. If --accessions is used, "
                    "each row should be the two accessions, separated by a tab. This will be used to nudge genomes that are haplotypes together")
    
    # ------------------------------------------------------------------
    # Snakemake execution options
    # ------------------------------------------------------------------
    smk = p.add_argument_group("snakemake options")
    smk.add_argument("--kmc", help="Run optional minimizer stats using KMC - requires KMC3 to be installed.",
                                 action="store_true")
    smk.add_argument("--threads", type=int, default=12,
                     help="Number of threads to use.")
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
        "genomes": str(Path(args.genomes).resolve()) if args.genomes else "",
        "prefix":  args.prefix,
        "name_conversions": args.name_conversions,
        "make_tree": args.make_tree,
        "mt_genomes": args.mt_genomes if args.mt_genomes else "",
        "mt_name_conversions": args.mt_name_conversions if args.mt_name_conversions else "",
        "date":              f"{args.prefix}_assemblies",
        "fpr":               args.fpr,
        "ntsynt_k":          args.ntsynt_k,
        "ntsynt_w":          args.ntsynt_w,
        "treefile":          args.tree,
        "ntsynt_viz_ribbon_adjust": args.ntsynt_viz_ribbon_adjust,
        "scripts_dir":     str(SCRIPTS_DIR.resolve()),
        "block_size": args.block_size if args.block_size else "",
        "indel":      args.indel if args.indel else "",
        "merge":      args.merge if args.merge else "",
        "scale": args.scale if args.scale else 100e6,
        "seq_length": args.seq_length if args.seq_length else "",
        "w_rounds": args.w_rounds if args.w_rounds else "",
        "hashes":        args.hashes if args.hashes is not None else "",
        "target_genome": args.target_genome,
        "viz_length":    args.viz_length if args.viz_length is not None else "",
        "viz_format":    args.format,
        "viz_width":     args.width if args.width is not None else "",
        "viz_dpi":       args.dpi if args.dpi is not None else "",
        "centromeres":   args.centromeres,
        "viz_order":     args.order,
        "viz_keep":      " ".join(args.keep) if args.keep else "",
        "no_arrow":      args.no_arrow,
        "optimize_ordering": args.optimize_ordering if args.optimize_ordering else "",
        "haplotypes":    args.haplotypes if args.haplotypes else "",
        "keep_incomplete":  args.keep_incomplete if args.keep_incomplete else "",
    }


def validate_paths(args: argparse.Namespace) -> None:
    """Abort early if required input files are missing."""
    errors = []
    if args.genomes and not Path(args.genomes).exists():
        errors.append(f" --genomes: file not found: {args.genomes}")
    elif args.genomes:
        with open(args.genomes, 'r', encoding="utf-8") as fin:
            for genome in fin:
                if not Path(genome.strip()).exists():
                    errors.append(f"Genome file listed in --genome not found: {genome}")
    if args.accessions and not Path(args.accessions).exists():
        errors.append(f" --accessions: file not found: {args.accessions}")
    if args.tree and not Path(args.tree).exists():
        errors.append(f"  --tree: file not found: {args.tree}")
    if args.make_tree and args.mt_genomes and not Path(args.mt_genomes).exists():
        errors.append(f"  --mt-genomes: file not found: {args.mt_genomes}")
    if args.make_tree and args.mt_genomes and not Path(args.mt_name_conversions).exists():
        errors.append(f"  --mt-name-conversions: file not found: {args.mt_name_conversions}")
    if args.centromeres and not Path(args.centromeres).exists():
        errors.append(f"  --centromeres: file not found: {args.centromeres}")
    if args.order and not Path(args.order).exists():
        errors.append(f"  --order: file not found: {args.order}")
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
        "--cores",     str(args.threads),
        "--printshellcmds",
        "--nolock",
    ]

    if args.kmc:
        cmd += ["all", "kmc"]

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
