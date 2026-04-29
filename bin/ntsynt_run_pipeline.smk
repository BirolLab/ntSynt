#!/usr/bin/env snakemake -s
import shutil
import sys
import os

# onsuccess:
#     shutil.rmtree(".snakemake", ignore_errors=True)

# Read in parameters
references = config["references"] if "references" in config else "Must specify 'references'"
k = config["kmer"] if "kmer" in config else 24
w = config["window"] if "window" in config else 1000
fpr = config["fpr"] if "fpr" in config else 0.025
hashes = config["hashes"] if "hashes" in config else 3
max_threads = config["threads"] if "threads" in config else 4
prefix = config["prefix"] if "prefix" in config else "ntSynt_out"
common = config["common"] if "common" in config else True
repeat = config["repeat"] if "repeat" in config else False
w_rounds = config["w_rounds"] if "w_rounds" in config else [100, 10]
indel_merge = config["indel_merge"] if "indel_merge" in config else 500
collinear_merge = config["collinear_merge"] if "collinear_merge" in config else "3w"
simplify_graph = config["simplify_graph"] if "simplify_graph" in config else True
benchmark = config["benchmark"] if "benchmark" in config else False
dev = config["dev"] if "dev" in config else False
min_block_size = config["block_size"] if "block_size" in config else 500

# If want to benchmark, use memusg or /usr/bin/time
# The snakemake benchmarking is quite inaccurate for shorter runtimes
benchmark_path = "" 
if benchmark:
    if benchmark_memusg := shutil.which("memusg"):
        benchmark_path = f"{benchmark_memusg} -t"
    elif benchmark_time := shutil.which("time"):
        benchmark_path = f"{benchmark_time} -v"
    else:
        print("WARNING: memusg and time not found in PATH. Benchmarks will not be tallied.")


# Get path to the base directory (where snakemake file and other scripts are located)
script_path = workflow.basedir

# basenames for input fasta files
config["dict_references"] = {os.path.basename(genome): genome for genome in references}

rule all:
    input: expand("{ref}.k{k}.w{w}.tsv", ref=config["dict_references"].keys(), k=k, w=w),
            f"{prefix}.synteny_blocks.tsv"

rule faidx:
    input: fa=lambda wildcards: config['dict_references'][wildcards.file]
    output: "{file}.fai"
    params: benchmarking=expand("{benchmark_path} -o {{file}}.faidx.time", benchmark_path=benchmark_path) if benchmark else []
    threads: 1
    shell: "{params.benchmarking} samtools faidx -o {output} {input.fa}"

rule make_common_bf:
    input: refs=references
    output: f"{prefix}.common.bf"
    threads: max_threads
    params: options=f"-p {prefix}.common --fpr {fpr} --hashes {hashes} -k {k}",
            path_to_script=f"{script_path}/ntsynt_make_common_bf",
            benchmarking=f"{benchmark_path} -o {prefix}.make_common_bf.time" if benchmark else []
    shell: "{params.benchmarking} {params.path_to_script} --genome {input.refs} {params.options} -t {threads}"

# For posterity, included but this is experimental
rule make_repeat_bf:
    input: refs=references
    output: expand("{prefix}.repeat.bf", prefix=prefix)
    threads: max_threads
    params: options=f"-p {prefix}.repeat --fpr {fpr} -k {k}",
            path_to_script=f"{script_path}/ntsynt_make_repeat_bfs.py",
            benchmarking=expand("{benchmark_path} -o {prefix}.make_repeat_bf.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} {params.path_to_script} --genome {input.refs} {params.options} -t {threads}"

rule indexlr:
    input: fa=lambda wildcards: config['dict_references'][wildcards.fasta],
            common=f"{prefix}.common.bf" if common is True else [],
            repeat=f"{prefix}.repeat.bf" if repeat is True else []
    output: expand("{{fasta}}.k{k}.w{w}.tsv", k=k, w=w)
    threads: 5
    resources: load=1
    params: options=f"-k {k} -w {w} --long --seq --pos",
            bf="-s" if common is True else [],
            repeat="-r" if repeat is True else [],
            benchmarking=expand("{benchmark_path} -o {{fasta}}.indexlr.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} indexlr {params.options} -t {threads} {params.bf} {input.common} {params.repeat} {input.repeat} {input.fa} > {output}"

rule ntsynt_synteny:
    input: mx=expand("{fasta}.k{k}.w{w}.tsv", fasta=config['dict_references'].keys(), k=k, w=w),
            common=f"{prefix}.common.bf" if common is True else [],
            repeat=f"{prefix}.repeat.bf" if repeat is True else [],
            fais=expand("{fasta}.fai", fasta=config['dict_references'].keys()),
            fastas=references
    output: f"{prefix}.synteny_blocks.tsv"
    threads: max_threads
    params: path_to_script=f"{script_path}/ntsynt_run.py",
            options=f"-k {k} -w {w} --w-rounds {w_rounds} -p {prefix} --bp {indel_merge} --collinear-merge {collinear_merge} -z {min_block_size}",
            common_bf="--common" if common is True else [],
            repeat_bf="--repeat" if repeat is True else [],
            simplify_graph="--simplify-graph" if simplify_graph is True else [],
            dev="--dev" if dev is True else [],
            benchmarking=f"{benchmark_path} -o {prefix}.synteny_blocks.time" if benchmark else [] 
    shell: "{params.benchmarking} python3 {params.path_to_script} {input.mx} {params.options} {params.common_bf} {input.common}  {params.simplify_graph} \
             --btllib_t {threads}  --fastas {input.fastas} {params.repeat_bf} {input.repeat} {params.dev}"
