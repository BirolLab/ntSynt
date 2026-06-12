#!/usr/bin/env snakemake -s
"""
phylogeny.smk  —  phylogenetic tree construction for ntSynt pipeline
=====================================================================
Constructs a phylogenetic tree from mitochondrial or nuclear sequences,
for use alongside the ntSynt-viz ribbon plot.

 1. If user supplies mitochondrial sequences (--mt-source user-fasta  +  --mt-fasta <path>),
    these sequences are used to construct a ML tree with mafft and iqtree
 2. If --mt-source auto is specified, the sequence reports for the downloaded accessions are 
    queried for embedded mitochondrial sequences. If found, continue with ML tree as with (1).
    Otherwise, construct tree using nuclear genomes
 3. If --mt-source nuclear or mt genomes cannot be found with --mt-source auto, use mashtree
    to construct a neighbour-joining tree directly from mash skecthes of the genomes.

Override flags (set in run_phylogeny.py driver):
  mt_source : "auto" | "user-fasta" | "nuclear"
  mt_fasta  : path to user-supplied MT FASTA  (required when mt_source=user-fasta)
"""

import os
import glob
import csv
from pathlib import Path
import re
import shutil

# ---------------------------------------------------------------------------
# Config aliases
# ---------------------------------------------------------------------------
GROUP       = config["taxonomic_group"]
FAM_LOW     = GROUP.lower()
DATE        = config["date"]
ASSEMBLY_DIR = f"{DATE}_assemblies"
NCBI_DATA_DIR = f"{ASSEMBLY_DIR}/{FAM_LOW}_assemblies/ncbi_dataset/data"

SEQ_REPORT    = config["seq_report"]        # from master pipeline
NAME_CONV     = config["name_conversion"]   # from master pipeline
FASTA_LIST    = config["fasta_list"]        # from master pipeline

MT_DIR        = f"phylogeny/mt"
TREE_DIR      = f"phylogeny/tree"

MT_SOURCE     = config.get("mt_source", "auto")   # auto | user-fasta | nuclear
MT_FASTA_USER = config.get("mt_fasta", "")        # only used when mt_source=user-fasta

SCRIPTS       = config.get("scripts_dir", {})
THREADS       = config.get("threads", 12)
PREFIX        = config.get("prefix", FAM_LOW)

# ---------------------------------------------------------------------------
# Determine mt_source at parse time when set to "auto":
# inspect the sequence report and check whether all species have MT.
# Sets MT_SOURCE_RESOLVED to "embedded", "download", or "nuclear".
# When mt_source != "auto" the user's choice is used directly.
# ---------------------------------------------------------------------------

def resolve_mt_source() -> str:
    """
    Inspect the sequence report to determine MT availability.
    Returns one of: "embedded", "download", "nuclear", "user-fasta".
    """
    if MT_SOURCE == "user-fasta":
        if not MT_FASTA_USER:
            raise ValueError("mt_source=user-fasta requires mt_fasta to be set.")
        return "user-fasta"

    if MT_SOURCE == "nuclear":
        return "nuclear"

    # MT_SOURCE == "auto": check the sequence report
    if not os.path.exists(SEQ_REPORT):
        # Report not yet generated — can't determine; default to nuclear
        raise ValueError(f"{SEQ_REPORT} does not exist - exiting.")

    # Load name conversion: fasta filename -> species name
    species_set = set()
    if os.path.exists(NAME_CONV):
        with open(NAME_CONV) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    fasta_base = parts[0]
                    # Strip suffix to get bare accession
                    m = re.match(r"(GC[AF]_\d+\.\d+)", fasta_base)
                    acc = m.group(1) if m else fasta_base
                    species_set.add(acc)

    total_assemblies = len(species_set)

    # Count assemblies with MT in the sequence report
    assemblies_with_mt = set()
    with open(SEQ_REPORT, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mol = row.get("Molecule type", row.get("mol-type", "")).strip()
            if mol.lower() in ("mitochondrion", "mitochondrial"):
                assemblies_with_mt.add(row.get("Assembly Accession", "").strip())
            mol = row.get("Chromosome name", "")
            if mol.lower() in ("mit", "mt"):
                assemblies_with_mt.add(row.get("Assembly Accession", "").strip())

    if len(assemblies_with_mt) == total_assemblies:
        print(f"  MT source: all {total_assemblies} assemblies have embedded MT.")
        return "embedded"
    elif len(assemblies_with_mt) > 0:
        print(
            f"  MT source: {len(assemblies_with_mt)}/{total_assemblies} assemblies have embedded MT."
        )
        return "download"
    else:
        print(f"  MT source: no embedded MT found.")
        return "download"


MT_SOURCE_RESOLVED = resolve_mt_source()

# ---------------------------------------------------------------------------
# Determine which tree pipeline to use
# ---------------------------------------------------------------------------

def tree_file(wildcards):
    ckpt = checkpoints.collect_mt_fastas.get()
    flag = ckpt.output.nuclear_flag

    with open(flag) as fh:
        use_nuclear = fh.read().strip() == "1"

    if use_nuclear:
        return f"{MT_DIR}/nuclear_mashtree.nwk"
    else:
        return f"{MT_DIR}/mt_assemblies.aln.fa.treefile"

# ---------------------------------------------------------------------------
# Final tree file (used by rule all)
# ---------------------------------------------------------------------------

FINAL_TREE = f"{TREE_DIR}/{PREFIX}.nwk"

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------

rule all:
    input:
        FINAL_TREE,


# ===========================================================================
# Look for mitochondrial sequences embedded in the sequence reports, 
# reverting to older assembly builds if needed.
# ===========================================================================

rule find_missing_mt_species:
    """
    Compare the name conversion file against the sequence report to find
    species that lack an MT sequence in their assembly.
    Writes a TSV of species names to query other assembly versions.
    """
    input:
        seq_report = SEQ_REPORT,
        name_conv  = NAME_CONV,
    output:
        missing_tsv = f"{MT_DIR}/download/missing_mt_species.txt",
    log:
        f"{MT_DIR}/download/find_missing.log",
    run:
        import csv, os

        # Assemblies that already have MT
        assemblies_with_mt = set()
        with open(input.seq_report, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                mol = row.get("Molecule type", row.get("mol-type", "")).strip().lower()
                if mol in ("mitochondrion", "mitochondrial"):
                    assemblies_with_mt.add(row.get("Assembly Accession", "").strip())
                mol = row.get("Chromosome name", "")
                if mol in ("MIT", "MT"):
                    assemblies_with_mt.add(row.get("Assembly Accession", "").strip())

        # All assemblies from name conversion: fasta_basename -> species_name
        missing_species = []
        with open(input.name_conv) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                fasta_base, species = parts[0], parts[1]
                # Extract accession from fasta basename (GCA_xxx.x prefix)
                m = re.match(r"(GC[AF]_\d+\.\d+)", fasta_base)
                acc = m.group(1) if m else fasta_base
                if acc not in assemblies_with_mt:
                    missing_species.append((acc, species))
                    print(f"  Missing MT for {acc} ({species})", file=open(log[0], "a"))

        os.makedirs(os.path.dirname(output.missing_tsv), exist_ok=True)
        with open(output.missing_tsv, "w") as fh:
            for acc, species in missing_species:
                fh.write(f"{acc}\t{species}\n")

        print(f"  {len(missing_species)} species missing MT sequences.")


rule download_missing_mt:
    """
    For each species missing an MT sequence, check all assembly versions for
    MT sequences via the sequence report. If found, download those assemblies
    and extract the MT sequence. If still not found, set the nuclear fallback flag.
    """
    input:
        missing_tsv = f"{MT_DIR}/download/missing_mt_species.txt",
    output:
        done          = f"{MT_DIR}/download/download_mt.done",
        fallback_flag = f"{MT_DIR}/download/use_nuclear.flag",
        missing_accs  = f"{MT_DIR}/download/missing_mt_accessions.txt",
    params:
        outdir    = f"{MT_DIR}/download",
        seq_report = f"{MT_DIR}/download/missing_seq-reports.tsv",
        mt_accs   = f"{MT_DIR}/download/missing_mt_found_accessions.txt",
    log:
        f"{MT_DIR}/download/download_mt.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        use_nuclear=0

        # If no species are missing MT, exit cleanly
        if [ ! -s {input.missing_tsv} ]; then
            echo "  No missing MT species — moving to next step." | tee -a {log}
            touch {output.missing_accs}
            echo "$use_nuclear" > {output.fallback_flag}
            touch {output.done}
            exit 0
        fi

        # ------------------------------------------------------------------
        # 1. Get the original accessions for missing species
        # ------------------------------------------------------------------
        cut -f1 {input.missing_tsv} > {output.missing_accs}

        # ------------------------------------------------------------------
        # 2. Query sequence reports across ALL assembly versions
        # ------------------------------------------------------------------
        echo "  Querying sequence reports for all assembly versions to attempt MT sequence retrieval..." | tee -a {log}
        datasets summary genome accession \
            --inputfile {output.missing_accs} \
            --report sequence \
            --assembly-level chromosome,complete \
            --assembly-version all \
            --as-json-lines \
        | dataformat tsv genome-seq \
            --fields accession,chr-name,genbank-seq-acc,mol-type,role \
        > {params.seq_report}

        # ------------------------------------------------------------------
        # 3. Check which original accessions have MT in any version
        # ------------------------------------------------------------------
        mlr --tsv \
            filter '${{Molecule type}} == "Mitochondrion" || ${{Molecule type}} == "Mitochondrial" || ${{Chromosome name}} == "MT" || ${{Chromosome name}} == "MIT";' \
            then cut -f "Assembly Accession" \
            then uniq -g "Assembly Accession" \
            {params.seq_report} \
        | tail -n +2 > {params.mt_accs}

        n_found=$(wc -l < {params.mt_accs})
        n_missing=$(wc -l < {output.missing_accs})
        echo "  Found MT in ${{n_found}}/${{n_missing}} missing assemblies (across all versions)." \
            | tee -a {log}

        # ------------------------------------------------------------------
        # 4. If not all species have MT, fall back to nuclear
        #    Otherwise, batch download + extract all at once
        # ------------------------------------------------------------------
        if [ "${{n_found}}" -ne "${{n_missing}}" ]; then
            echo "  WARNING: MT not found for all missing species — falling back to nuclear." \
                | tee -a {log}
            # Log which accessions are missing
            comm -23 \
                <(sort {output.missing_accs}) \
                <(sort {params.mt_accs}) \
            | while read -r acc; do
                echo "    No MT found for: $acc" | tee -a {log}
            done
            use_nuclear=1
        else
            echo "  MT found for all missing species — downloading in batch..." | tee -a {log}

            datasets download genome accession \
                --inputfile {params.mt_accs} \
                --include genome \
                --assembly-level complete,chromosome \
                --filename {params.outdir}/missing_mt_assemblies.zip \
                2>> {log}

            unzip -o {params.outdir}/missing_mt_assemblies.zip \
                  -d {params.outdir}/missing_mt_assemblies 2>> {log}
            rm -f {params.outdir}/missing_mt_assemblies.zip

            # Extract MT sequence for each accession into <acc>.mt.fa
            while IFS=$'\t' read -r acc; do
                mt_gbk=$(mlr --tsv \
                             filter -s my_acc="$acc" \
                                 '${{Assembly Accession}} == @my_acc && \
                                 (${{Molecule type}} == "Mitochondrion" || ${{Molecule type}} == "Mitochondrial" || ${{Chromosome name}} == "MT" || ${{Chromosome name}} == "MIT")' \
                             then cut -f "GenBank seq accession" \
                             {params.seq_report} \
                         | tail -n +2 | head -n 1)

                fna=$(find {params.outdir}/missing_mt_assemblies \
                           -path "*/${{acc}}/*.fna" | head -n1)

                if [ -z "$fna" ]; then
                    echo "  WARNING: no .fna found for $acc after download." | tee -a {log}
                    use_nuclear=1
                    continue
                fi

                out={params.outdir}/${{acc}}.mt.fa
                seqtk subseq $fna <(echo $mt_gbk) > "$out"
                echo "  Extracted MT $mt_gbk for $acc -> $out" | tee -a {log}
            done < {params.mt_accs}
        fi

        echo "$use_nuclear" > {output.fallback_flag}
        touch {output.done}
        """
# ===========================================================================
# Extract MT sequences from chromosome assemblies
# ===========================================================================

rule extract_mt_from_assemblies:
    """
    Extract mitochondrial sequences from the downloaded chromosome FASTAs
    using the sequence report (Molecule type == Mitochondrion | Mitochondrial OR Chromosome name == MT | MIT).
    Writes one <accession>.mt.fa per assembly into MT_DIR/embedded/.
    """
    input:
        seq_report = SEQ_REPORT,
        fasta_list = FASTA_LIST,
        missing_done = rules.download_missing_mt.output.missing_accs,
        nuclear_fallback = rules.download_missing_mt.output.fallback_flag,
    output:
        done = f"{MT_DIR}/embedded/extract_mt.done",
    params:
        outdir     = f"{MT_DIR}/embedded",
        fasta_root = NCBI_DATA_DIR,
    log:
        f"{MT_DIR}/embedded/extract_mt.log",
    shell:
        r"""
        set -euxo pipefail
        mkdir -p {params.outdir}

        # Check nuclear fallback flag
        nuclear_flag=$(cat {input.nuclear_fallback})
        if [ "$nuclear_flag" == "1" ]; then
            # Don't extract mt assemblies if going to do nuclear mashtree anyway
            echo "Skipping extracting mt from assemblies, as nuclear flag set."
            touch {output.done}
            exit 0
        fi

        # Iterate over assemblies that have MT in the sequence report
        mlr --tsv \
            filter '$["Molecule type"] == "Mitochondrion" || $["Molecule type"] == "Mitochondrial" || $["Chromosome name"] == "MT" || $["Chromosome name"] == "MIT"' \
            then cut -f "Assembly Accession,GenBank seq accession" \
            {input.seq_report} \
        | tail -n +2 \
        | while IFS=$'\t' read -r acc gbk_acc; do
            fasta=$(find {params.fasta_root}/${{acc}} -name "*.fna" 2>/dev/null | head -n1)
            if [ -z "$fasta" ]; then
                echo "  WARNING: no .fna found for $acc, skipping." | tee -a {log}
                continue
            fi
            out={params.outdir}/${{acc}}.mt.fa
            echo "  Extracting MT ${{gbk_acc}} from $acc -> $out" | tee -a {log}
            seqtk subseq $fasta  <(echo $gbk_acc) > $out
        done

        touch {output.done}
        """


# ===========================================================================
# Collect all extracted/downloaded MT FASTAs (or user-supplied)
# and decide whether to proceed with MT or fall back to nuclear
# ===========================================================================

checkpoint collect_mt_fastas:
    """
    Gather all per-assembly MT FASTAs (from embedded extraction and/or
    download) into a single directory and write a file-of-files.
    Checks the fallback flag — if nuclear is needed, writes a separate flag.
    """
    input:
        embedded_done  = f"{MT_DIR}/embedded/extract_mt.done"
                         if MT_SOURCE_RESOLVED in ("embedded", "download") else [],
        download_done  = f"{MT_DIR}/download/download_mt.done"
                         if MT_SOURCE_RESOLVED == "download" else [],
        fallback_flag  = f"{MT_DIR}/download/use_nuclear.flag"
                         if MT_SOURCE_RESOLVED == "download" else [],
    output:
        fof          = f"{MT_DIR}/mt_assemblies.fof",
        nuclear_flag = f"{MT_DIR}/use_nuclear.flag",
    params:
        embedded_dir = f"{MT_DIR}/embedded",
        download_dir = f"{MT_DIR}/download",
        mt_fasta     = MT_FASTA_USER,
        mt_source    = MT_SOURCE_RESOLVED,
    log:
        f"{MT_DIR}/collect_mt.log",
    run:
        use_nuclear = False
        mt_files = []

        if params.mt_source == "user-fasta":
            mt_files = [params.mt_fasta]
            print(f"  Using user-supplied MT FASTA: {params.mt_fasta}")

        elif params.mt_source == "nuclear":
            use_nuclear = True

        else:
            # Collect from embedded and/or download directories
            mt_files += sorted(glob.glob(f"{params.embedded_dir}/*.mt.fa"))
            if params.mt_source == "download":
                mt_files += sorted(glob.glob(f"{params.download_dir}/*.mt.fa"))

                # Check fallback flag from download rule
                with open(input.fallback_flag) as fh:
                    if fh.read().strip() == "1":
                        print(
                            "  WARNING: not all species have MT sequences — "
                            "falling back to nuclear genome tree.",
                            file=open(log[0], "a")
                        )
                        use_nuclear = True

        os.makedirs(os.path.dirname(output.fof), exist_ok=True)

        with open(output.fof, "w") as fh:
            for f in mt_files:
                fh.write(os.path.abspath(f) + "\n")

        with open(output.nuclear_flag, "w") as fh:
            fh.write("1" if use_nuclear else "0")

        print(f"  Collected {len(mt_files)} MT FASTAs. Use nuclear: {use_nuclear}")


# ===========================================================================
# MT tree: mafft + iqtree
# ===========================================================================

rule concat_mt_fastas:
    """Concatenate all per-assembly MT FASTAs into one multi-FASTA."""
    input:
        fof = f"{MT_DIR}/mt_assemblies.fof",
    output:
        fa  = f"{MT_DIR}/mt_assemblies.fa",
    shell:
        r"""
        cat $(cat {input.fof} | tr "\n" " ") > {output.fa}
        """


rule mafft_align:
    input:
        fa  = f"{MT_DIR}/mt_assemblies.fa",
    output:
        aln = f"{MT_DIR}/mt_assemblies.aln.fa",
    params:
        maxiterate = 1000,
    threads: THREADS
    log:
        f"{MT_DIR}/mafft.log",
    shell:
        r"""
        mafft --maxiterate {params.maxiterate} --thread {threads} \
            {input.fa} > {output.aln} 2> {log}
        """


rule iqtree:
    input:
        aln = f"{MT_DIR}/mt_assemblies.aln.fa",
    output:
        tree = f"{MT_DIR}/mt_assemblies.aln.fa.treefile",
    threads: THREADS
    log:
        f"{MT_DIR}/iqtree.log",
    shell:
        r"""
        iqtree -s {input.aln} -nt {threads} &> {log}
        """


# ===========================================================================
# Nuclear tree: mashtree
# ===========================================================================

rule mashtree:
    input:
        fasta_list = FASTA_LIST,
    output:
        tree = f"{MT_DIR}/nuclear_mashtree.nwk",
    threads: THREADS
    shell:
        r"""
        set -eux -o pipefail
        genome_size=$(cat $(head -n1 {input.fasta_list}) | awk '/^>/{{next}} {{size += length($0)}} END{{print size}}')

        # Build file-of-files from fasta list
        fastas=$(cat {input.fasta_list}  |tr "\n" " ")
        mashtree.pl --numcpus {threads} \
            --genomesize ${{genome_size}} --tempdir tmp \
            ${{fastas}} \
            > {output.tree} 
        """


# ===========================================================================
# Select the right treefile (MT or nuclear) and rename leaves
# ===========================================================================

rule make_mt_name_conversion:
    """
    Build a sequence-ID -> species name translation file for rename_newick.
    The treefile leaf labels will be MT sequence accessions (e.g. CM012345.1),
    not assembly accessions. This rule joins:
      seq_report:    GenBank seq accession -> Assembly Accession
      name_conv:     Assembly Accession (from filename) -> Species name
    to produce:
      GenBank seq accession -> Species name

    If the nuclear flag is set, leaf labels are assembly accessions (filename
    stems), so the output maps those directly to species names instead.
    """
    input:
        seq_report   = SEQ_REPORT,
        name_conv    = NAME_CONV,
        nuclear_flag = f"{MT_DIR}/use_nuclear.flag",
    output:
        mt_name_conv = f"{TREE_DIR}/mt_name_conversion.tsv",
    run:
        nuclear = open(input.nuclear_flag).read().strip() == "1"

        # Build accession -> species name from name conversion file
        # Filenames like GCA_965165685.3.chr.fa -> extract accession prefix
        acc_to_species = {}
        with open(input.name_conv) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                fasta_base, species = parts[0], parts[1]
                m = re.match(r"(GC[AF]_\d+\.\d+)", fasta_base)
                acc = m.group(1) if m else fasta_base
                acc_to_species[acc] = species

        os.makedirs(TREE_DIR, exist_ok=True)

        if nuclear:
            # Leaf labels are filename stems (assembly accessions) — map directly
            with open(output.mt_name_conv, "w") as out:
                for acc, species in acc_to_species.items():
                    out.write(f"{acc}.chr\t{species}\n")
        else:
            # Join via sequence report: GenBank seq accession -> assembly accession
            with open(input.seq_report, newline="") as fh, \
                 open(output.mt_name_conv, "w") as out:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    mol = row.get("Molecule type", row.get("mol-type", "")).strip().lower()
                    if mol not in ("mitochondrion", "mitochondrial"):
                        continue
                    gbk_acc = row.get("GenBank seq accession", "").strip()
                    asm_acc = row.get("Assembly Accession", "").strip()
                    species = acc_to_species.get(asm_acc, asm_acc)
                    if gbk_acc:
                        out.write(f"{gbk_acc}\t{species}\n")

rule rename_newick:
    input:
        tree=tree_file,
        name_conv = f"{TREE_DIR}/mt_name_conversion.tsv",
    output:
        renamed = FINAL_TREE,
    params: 
        script = f"{SCRIPTS}/rename_newick.py",
    log:
        f"{TREE_DIR}/rename_newick.log",
    shell:
        r"""
        python3 {params.script} {input.tree} {input.name_conv} \
            > {output.renamed} 2> {log}
        """


