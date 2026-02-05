#!/usr/bin/env python3
'''
Run the minimizer graph stage of ntSynt
'''
import argparse
from ntsynt_synteny import NtSyntSynteny

NTSYNT_VERSION = 'ntSynt v1.0.5'

def parse_arguments():
    "Parse arguments from argparse"
    parser = argparse.ArgumentParser(description="Run the dynamic minimizer graph stage of ntSynt")
    parser.add_argument("FILES", nargs="+", help="Minimizer TSV files of input assemblies")
    parser.add_argument("--fastas", nargs="+", help="Assembly fasta files", required=True, type=str)
    parser.add_argument("-n", help="Minimum edge weight [Number of input assemblies]", default=0, type=int)
    parser.add_argument("-p", help="Output prefix [out]",
                                default="out", type=str, required=False)
    parser.add_argument("-k", help="k-mer size used for minimizer step", required=True, type=int)
    parser.add_argument("-w", help="Window size used for minimizers", required=True, type=int)
    parser.add_argument("-z", help="Minimum synteny block size (bp) [500]", type=int, default=500)
    parser.add_argument("--filter", help="Type of repeat filtering", choices=["Filter", "Indexlr"], type=str)
    parser.add_argument("--common", help="Input common BF for minimizer selection", type=str)
    parser.add_argument("--repeat", help="Repeat BF (must be included if --filter is specified)", type=str)
    parser.add_argument("--btllib_t", help="Number of threads for btllib wrapper functions "\
                                "(computing minimizers, reading fasta file) [4]", type=int, default=4)
    parser.add_argument("--w-rounds", help="decreasing list of 'w' values to use for refining ends",
                                default=[100, 10], nargs="+", type=int)
    parser.add_argument("--bp", help="Maximum tolerated indel size [500]",
                                default=500, type=int)
    parser.add_argument("--collinear-merge", help="Maximum distance between collinear blocks for merging"\
                                "(length in bp or string in the form '<num>w' to indicate multiples of w) [1w]",
                                default="1w", type=str, required=False)
    parser.add_argument("--simplify-graph", help="Run minimizer graph simplification",
                                action="store_true")
    parser.add_argument('-m', help="Require at least m %% of minimizer positions to be "\
                                "increasing/decreasing to assign contig orientation [90]\n ",
                                default=90, type=int)
    parser.add_argument("--dev", action="store_true", help="Developer mode - retain intermediate files"\
                                                            "and more verbose logging")
    parser.add_argument("--interarrivals", action="store_true",
                                help="Output interarrival distances in initial graph")
    parser.add_argument("-v", "--version", action='version', version=NTSYNT_VERSION)

    return parser.parse_args()

def main():
    "Run ntSynt"
    print(f"Running {NTSYNT_VERSION}", flush=True)
    args = parse_arguments()
    NtSyntSynteny(args).main_synteny()

if __name__ == "__main__":
    main()
