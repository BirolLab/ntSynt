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
import os
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

SNAKEFILE = Path(__file__).parent / "Snakefile"
SCRIPTS_DIR = Path(__file__).parent / "scripts"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run the multi-genome synteny analysis pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ------------------------------------------------------------------
    # Required inputs
    # ------------------------------------------------------------------
    req = p.add_argument_group("required inputs")
    req.add_argument(
        "--accessions", required=True, metavar="TSV",
        help=(
            "TSV listing genome assemblies. Must contain an 'Assembly Accession' "
            "column. Optional 'Species Name' or 'Organism Name' column is used for "
            "ntSynt-viz display labels."
        ),
    )
    req.add_argument(
        "--group", required=True, metavar="NAME",
        help="Taxonomic group name, used as a prefix throughout (e.g. 'ichneumonidae').",
    )
    req.add_argument(
        "--tax-level", default="family", metavar="COLUMN",
        help=(
            "Column name in the accessions TSV to filter on "
            "(e.g. 'family', 'order', 'class'). Default: 'family'."
        ),
    )
    req.add_argument(
        "--tax-value", default=None, metavar="VALUE",
        help=(
            "Value to match in --tax-level column (case-sensitive). "
            "Defaults to the capitalised form of --group if not supplied."
        ),
    )
    req.add_argument(
        "--block-stats-script", required=True, metavar="PATH",
        help="Path to denovo_synteny_block_stats.py (from the ntSynt repo).",
    )
    req.add_argument(
        "--mx-stats-script", required=True, metavar="PATH",
        help="Path to analyze_mx.py (minimizer stats script).",
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
    opt.add_argument("--tree", default="", metavar="NEWICK",
                     help="Optional Newick tree file for ntSynt-viz. Omit to skip.")
    opt.add_argument(
        "--date", required=True, metavar="YYMMDD",
        help="Date string used to name the download directory (e.g. '250528').",
        default=datetime.datetime.now().strftime("%y%m%d"),
        required=False
    )
    opt.add_argument("--ntsynt-viz_ribbon-adjust", type=float, default=0.2,
                     help="Adjustment factor for ntSynt-viz ribbons. Increase if ribbon plot labels are cut off. [0.2]")
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

    return p.parse_args()


def build_config(args: argparse.Namespace) -> dict:
    """Translate parsed CLI args into the config dict the Snakefile expects."""
    return {
        "accessions_tsv":    str(Path(args.accessions).resolve()),
        "taxonomic_group":   args.group,
        "tax_level":         args.tax_level,
        "tax_value":         args.tax_value if args.tax_value else args.group.capitalize(),
        "date":              args.date,
        "fpr":               args.fpr,
        "ntsynt_k":          args.ntsynt_k,
        "ntsynt_w":          args.ntsynt_w,
        "treefile":          args.tree,
        "ntsynt_viz_ribbon_adjust": args.ntsynt_viz_ribbon_adjust,
        "scripts": {
            "block_stats":              str(Path(args.block_stats_script).resolve()),
            "mx_stats":                 str(Path(args.mx_stats_script).resolve()),
            # Bundled scripts — always resolved relative to this driver file
            "extract_synteny_regions":  str(SCRIPTS_DIR / "extract_synteny_regions.py"),
            "plot_dists_r":             str(SCRIPTS_DIR / "plot_dists.R"),
        },
    }


def validate_paths(args: argparse.Namespace) -> None:
    """Abort early if required input files are missing."""
    errors = []
    for label, path in [
        ("--accessions",        args.accessions),
        ("--block-stats-script", args.block_stats_script),
        ("--mx-stats-script",   args.mx_stats_script),
    ]:
        if not Path(path).exists():
            errors.append(f"  {label}: file not found: {path}")
    if args.tree and not Path(args.tree).exists():
        errors.append(f"  --tree: file not found: {args.tree}")
    if errors:
        print("ERROR: the following required files were not found:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    validate_paths(args)
    config = build_config(args)

    # Import snakemake here so the script is importable even if snakemake is
    # not installed (e.g. during testing / --help).
    try:
        import snakemake as smk_module
    except ImportError:
        print("ERROR: snakemake Python package is not importable in this environment.",
              file=sys.stderr)
        sys.exit(1)

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

    success = smk_module.snakemake(
        snakefile   = str(SNAKEFILE),
        config      = config,
        cores       = args.cores,
        dryrun      = args.dry_run,
        forcerun    = args.forcerun or [],
        until       = args.until or [],
        printshellcmds = True,
        # Surface snakemake's own output directly to the terminal
        quiet       = False,
    )

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
