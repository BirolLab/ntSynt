# ntSynt Multi-Genome Synteny Pipeline

A Snakemake pipeline for end-to-end multi-genome synteny analysis using
[ntSynt](https://github.com/bcgsc/ntSynt). Starting from a TSV of genome
accessions, the pipeline downloads assemblies, measures pairwise divergence,
runs synteny block detection, and produces a suite of summary statistics and
visualisations.

---

## Table of contents

1. [Dependencies](#dependencies)
2. [Input TSV format](#input-tsv-format)
3. [Pipeline structure](#pipeline-structure)
4. [Usage](#usage)
5. [Outputs](#outputs)
6. [Notes](#notes)

---

## Dependencies

* Python 3.10+
* Snakemake
* [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
* samtools
* bedtools
* mash
* abyss
* miller
* ntSynt
* ntSynt-viz
* seqtk
* weasyprint
* mafft
* mashtree 1.4.6 + libnsl
* iqtree
* entrez-direct

---

## Input TSV format

The input TSV is generated using the NCBI `datasets` CLI. It must contain one
assembly per row and include the following columns (additional columns are
ignored):

### Required columns

| Column | Description |
|--------|-------------|
| `Assembly Accession` | NCBI accession (e.g. `GCA_964059415.1`) |
| `Assembly Level` | e.g. `Complete Genome`, `Chromosome`, `Scaffold` |
| Taxonomic rank columns | At minimum those used for filtering, e.g. `family`, `order`, `class` |

### Optional columns (used for display labels in ntSynt-viz)

| Column | Description |
|--------|-------------|
| `Species Name` | Used as the display name if present |
| `Organism Name` | Fallback if `Species Name` is absent |

### Generating the TSV with NCBI datasets

```bash
# 1. Search for assemblies and download a summary
datasets summary genome taxon "Lucinidae" \
    --assembly-level complete,chromosome \
    --as-json-lines \
    | dataformat tsv genome \
        --fields accession,assminfo-level,organism-name,organism-tax-id \
    > lucinidae_assemblies.tsv

# 2. Add taxonomic rank columns (order, family, genus, etc.) if not already present.
#    These can be joined from the NCBI taxonomy database or a pre-built ranks TSV.
```

### Pre-filtering recommendations

Before running the pipeline, filter the TSV to:

- **One assembly per species** — the pipeline does not deduplicate automatically.
  Retaining multiple assemblies per species will inflate divergence estimates and
  ntSynt runtime. Keep the highest-quality assembly (prefer `Complete Genome` over
  `Chromosome`, and higher contig N50).

Example filtering with miller:

```bash
mlr --tsv \
    then uniq -f "species" \
    assemblies_raw.tsv > assemblies_filtered.tsv
```

---

## Usage

The pipeline is launched via `run_pipeline.py`. 

```
python ntsynt_pipeline.py --help
```

### Minimal example

```bash
python run_pipeline.py \
    --accessions lucinidae_assemblies.tsv \
    --group      lucinidae \
    --tax-level  family \
    --tax-value  Lucinidae \
```


### Full argument reference

| Argument | Default | Description |
|----------|---------|-------------|
| `--accessions` | — | Input assemblies TSV (required) |
| `--group` | — | Taxonomic group name, used as output prefix (required) |
| `--tax-level` | `family` | TSV column name to filter on |
| `--tax-value` | capitalised `--group` | Value to match in `--tax-level` column (case-insensitive) |
| `--date` | today (`YY-MM-DD`) | Date string used to name the download directory |
| `--fpr` | `0.025` | Bloom filter false-positive rate for ntSynt |
| `--ntsynt-k` | `24` | k-mer size for ntSynt |
| `--ntsynt-w` | `1000` | Minimizer window size for ntSynt |
| `--tree` | _(none)_ | Newick tree file for ntSynt-viz (optional) |
| `--cores` | `12` | CPU cores available to Snakemake |
| `--dry-run` / `-n` | — | Preview rules without executing |
| `--forcerun` | — | Force re-execution of one or more named rules |
| `--until` | — | Stop pipeline after the named rule(s) |
| `--snakemake-args` | — | Additional flags passed verbatim to Snakemake |

---

## Outputs

All outputs are written relative to the working directory from which the
pipeline is launched.

```
<date>_assemblies/
├── <group>_biogenomes.tsv             # filtered input TSV (this taxonomic group only)
├── <group>_accessions.txt             # accession list passed to datasets
├── <group>_fasta_list.txt             # flat list of all downloaded FASTA paths
├── <group>_name_conversion.tsv        # filename → display name map for ntSynt-viz
├── <group>_abyss_fac.tsv              # chromosome/contig length stats
└── <group>_assemblies/
    └── ncbi_dataset/data/<accession>/*.fna   # downloaded genome FASTAs + .fai indices

mash/
├── full_assemblies_dist.tsv           # all-vs-all pairwise mash distances (full genomes)
├── max_mash_dist.txt                  # scalar max distance used as ntSynt -d
├── syntenic_dists.tsv                 # pairwise distances across syntenic regions
├── non_syntenic_dists.tsv             # pairwise distances across non-syntenic regions
├── mash_divergence_boxplot.pdf        # divergence distribution plot (all three categories)
└── mash_divergence_boxplot.png

ntsynt_run/
├── ntSynt.k24.w1000.synteny_blocks.tsv        # main ntSynt output
├── ntSynt.k24.w1000.synteny_blocks.stats.tsv  # block-level summary statistics
├── ntSynt.k24.w1000.mx.stats.tsv              # minimizer statistics
├── ntSynt.k24.w1000.discontinuity_reasons.tsv # tally of block break reasons
├── regions/
│   ├── <accession>.syntenic.fa        # syntenic sequence per assembly
│   └── <accession>.non_syntenic.fa    # non-syntenic sequence per assembly
└── ntsynt-viz/
    └── <group>.*                      # ntSynt-viz ribbon plots and supporting files
```

---

## Notes

**ntSynt divergence parameter** — the `-d` value passed to ntSynt is derived
automatically as the maximum pairwise mash distance across all full assemblies,
multiplied by 100 (to convert to the percentage scale ntSynt expects). No cap
is applied; if your group is very divergent you may want to inspect
`mash/max_mash_dist.txt` and override manually via `--snakemake-args`.

**Nested Snakemake / lock issues** — ntSynt itself runs an internal Snakemake
process. If you encounter a `LockException`, pass `--nolock` to the inner call
via the ntSynt `--nolock` flag (if supported), or ensure the ntSynt working
directory is set to a subdirectory distinct from the outer pipeline's working
directory.

**ntSynt-viz tree** — pass a Newick tree file with `--tree` to include a
phylogenetic guide tree in the ribbon plot. If omitted, ntSynt-viz runs without
a tree and assemblies are ordered by their appearance in the blocks TSV.
