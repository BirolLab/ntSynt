![GitHub Release](https://img.shields.io/github/v/release/BirolLab/ntsynt)
![Conda Downloads](https://img.shields.io/conda/dn/bioconda/ntsynt?label=Conda%20downloads)
[![Preprint](https://img.shields.io/badge/Preprint-blue.svg)](https://doi.org/10.1101/2024.02.07.579356)

![Logo](https://github.com/BirolLab/ntSynt/blob/main/logo/ntsynt-logo.png)

# ntSynt

Multi-genome macrosynteny detection using a dynamic minimizer graph approach.

## Contents
1. [Description of ntSynt](#description-of-ntsynt)
2. [Credits](#credits)
3. [Citing ntSynt](#citing-ntsynt)
4. [Usage](#usage)
5. [Installation instructions](#installation)
6. [Example](#example-command)
7. [Output files](#output-files)
8. [Example benchmarks](#example-benchmarks)
9. [Assessment](#basic-assessment-of-synteny-blocks)
10. [Tips / Visualization](#tips)
11. [License](#license)

## Description of ntSynt

ntSynt takes multiple genomes as input, and will compute synteny blocks that are in common with each of these input assemblies. These macrosyntenic blocks can enable a wide variety of comparative genomics studies between multiple genomes of varying divergences. ntSynt builds on the [ntJoin](https://github.com/BirolLab/ntJoin) codebase.

For more technical information about the various steps in ntSynt, see our [wiki page](https://github.com/BirolLab/ntSynt/wiki/Description-of-the-ntSynt-algorithm).

## Credits

Concept: Lauren Coombe and Rene Warren

Design and implementation: Lauren Coombe

## Citing ntSynt
If you use ntSynt in your work, please cite:

Lauren Coombe, Parham Kazemi, Johnathan Wong, Inanc Birol, René L. Warren. 2025. ntSynt: multi-genome synteny detection using minimizer graph mappings. BMC Biology. 23:367. https://doi.org/10.1186/s12915-025-02455-w

## Usage

```
usage: ntSynt [-h] [--fastas_list FASTAS_LIST] -d DIVERGENCE [-p PREFIX] [-k K] [-w W] [-t T] [--fpr FPR] [-b BLOCK_SIZE] [--merge MERGE]
              [--w_rounds W_ROUNDS [W_ROUNDS ...]] [--indel INDEL] [-n] [--benchmark] [-f] [--dev] [-v]
              [fastas ...]

ntSynt: Multi-genome synteny detection using minimizer graphs

positional arguments:
  fastas                Input genome fasta files

optional arguments:
  -h, --help            show this help message and exit
  --fastas_list FASTAS_LIST
                        File listing input genome fasta files, one per line
  -d DIVERGENCE, --divergence DIVERGENCE
                        Approx. maximum percent sequence divergence between input genomes (Ex. -d 1 for 1% divergence).
                        This will be used to set --indel, --merge, --w_rounds, --block_size
                        See below for default values - You can also set any of those parameters yourself, which will override these settings.
  -p PREFIX, --prefix PREFIX
                        Prefix for ntSynt output files [ntSynt.k<k>.w<w>]
  -k K                  Minimizer k-mer size [24]
  -w W                  Minimizer window size [1000]
  -t T                  Number of threads [12]
  --fpr FPR             False positive rate for Bloom filter creation [0.025]
  -b BLOCK_SIZE, --block_size BLOCK_SIZE
                        Minimum synteny block size (bp)
  --merge MERGE         Maximum distance between collinear synteny blocks for merging (bp). 
                        Can also specify a multiple of the window size (ex. 3w)
  --w_rounds W_ROUNDS [W_ROUNDS ...]
                        List of decreasing window sizes for synteny block refinement
  --indel INDEL         Threshold for indel detection (bp)
  -n, --dry-run         Print out the commands that will be executed
  --benchmark           Store benchmarks for each step of the ntSynt pipeline
  -f, --force           Run all ntSynt steps, regardless of existing output files
  --dev                 Run in developer mode to retain intermediate files, log verbose output
  -v, --version         show program's version number and exit
```
Given the approximate maximum divergence between the supplied genomes, ntSynt will set these default parameters:
|Divergence range|Default parameters|
|----|----|
|< 1%|--block_size 500 --indel 10000 --merge 10000 --w_rounds 100 10|
|1% - 10%|--block_size 1000 --indel 50000 --merge 100000 --w_rounds 250 100|
|>10%|	--block_size 10000 --indel 100000 --merge 1000000 --w_rounds 500 250|

Any of these parameters can be overridden by specifying them in your command. While these settings work generally well for the associated divergence range, we highly recommend customizing them for your particular requirements.

If not already known, we suggest using [Mash](https://github.com/marbl/Mash) to approximate the divergences, or mutation rates, between your compared genomes. Here are some example pairwise divergences as determined using Mash and reference genomes:
|Genomes|Approximate Mash sequence divergence (%)|
|---|---|
|human and chimpanzee|1.3|
|human and bonobo|1.3|
|chimapanzee and bonobo|0.5|
|mouse and rat|12|

## Installation

### Installing via conda
```
conda install -c bioconda -c conda-forge ntsynt
```

### Dependencies
- python 3.9+ with modules:
  - [intervaltree](https://github.com/chaimleib/intervaltree)
  - [pybedtools](https://daler.github.io/pybedtools/)
  - [ncls](https://github.com/pyranges/ncls)
  - [python-igraph](https://python.igraph.org/en/stable/)
- GCC 6+ or Clang 5+ (with OpenMP and C++17 support)
- [btllib v1.6.2+](https://github.com/bcgsc/btllib)
- [meson](https://mesonbuild.com/)
- [ninja](https://ninja-build.org/)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [samtools](http://www.htslib.org/)
- [seqtk](https://github.com/lh3/seqtk)

### Installing ntSynt from the source code
```
meson setup build --prefix=/path/to/desired/install/location
cd build
ninja install
```

### Testing ntSynt installation
Test your ntSynt installation using our provided demo:
```
cd tests
./run_ntSynt_demo.sh 
```
Once the script has executed successfully, you can compare the output files with those in tests/expected_results

### Example command
To compute the synteny blocks between 3 assemblies (assembly1.fa, assembly2.fa, assembly3.fa) with default parameters, where the maximum sequence divergence among these is ~5%, run:
```
ntSynt -d 5 assembly1.fa assembly2.fa assembly3.fa
```

### Output files
The main output file has the naming scheme `<prefix>.synteny_blocks.tsv`. This contains the synteny blocks computed in a TSV format.

The columns of this output synteny blocks TSV:
1. Synteny block ID - Lines with the same ID are part of the same synteny block
2. Genome file name
3. Genome chromosome/contig
4. Genome start coordinate
5. Genome end coordinate
6. Chromosome/contig strand
7. Number of mapped minimizers in this synteny block
8. Reason for discontinuity with previous synteny block

### Example benchmarks

|Synteny comparison|Number of genomes|Max. genome size (Gbp)|Wall clock time (min)|Peak memory usage (GB)|
|-----|-----|-----|-----|-----|
|Human genomes (0.1% simulated divergence)|2|3|26|34|
|Great ape genomes|4|3|48|32|
|Bee genomes|11|0.44|15|4|

### Basic assessment of synteny blocks
For a basic statistical summary of the computed synteny blocks, you can use the script `denovo_synteny_block_stats.py` found in `analysis_scripts`:
```
python3 denovo_synteny_block_stats.py -h
usage: denovo_synteny_block_stats.py [-h] --tsv TSV --fai FAI [FAI ...]

Compute de novo stats on synteny blocks

optional arguments:
  -h, --help           show this help message and exit
  --tsv TSV            ntSynt synteny block file
  --fai FAI [FAI ...]  FAI files for the compared genomes
```
More information can be found on our [wiki page](https://github.com/BirolLab/ntSynt/wiki/de-novo-statistics-summary)

### Tips / Visualization <a name=tips></a>
- To lower the peak memory usage, increase the false positive rate (--fpr) for the constructed Bloom filter
- Customize parameters such as --merge, --indel, --block_size and --w_rounds for your particular input data and research questions
- For visualizing the multi-genome output synteny blocks, please refer to 1) [ntSynt-viz](https://github.com/BirolLab/ntSynt-viz) and/or 2)the sub-directory [visualization_scripts](https://github.com/BirolLab/ntSynt/tree/main/visualization_scripts)
- If you do not know the approximate sequence divergence between the input assemblies, we recommend using [Mash](https://github.com/marbl/Mash) to estimate the divergences


## License
ntSynt Copyright (c) 2023-present British Columbia Cancer Agency Branch. All rights reserved.

ntSynt is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein prebstein@bccancer.bc.ca
