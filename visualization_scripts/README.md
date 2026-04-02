## NOTE: These example visualization scripts are kept here for posterity, but check out our visualization tool [ntSynt-viz](https://github.com/BirolLab/ntSynt-viz) for a more user-friendly and streamlined approach to generating synteny visualizations!

# Example visualization scripts for ntSynt

Here, we provide basic examples scripts for generating ribbon plots and chromosome sequence painting plots to visualize synteny blocks computed by ntSynt. Each R script script can be seen as a starting point, and can be customized as needed. 

## Ribbon plots
The R package gggenomes (https://thackl.github.io/gggenomes/) is used to generate ribbon plots to visualize multi-genome synteny blocks.

### Steps:
1. Format the ntSynt synteny blocks for the ribbon plots R script using the provided script. Before running `format_blocks_gggenomes.sh`, add the `visualization_scripts` directory to your PATH.
```
Usage: format_blocks_gggenomes.sh <synteny blocks TSV> <prefix> <length threshold> <assembly to use for colour> <FAI> <FAI> [FAI..]
```
The order of assembly FAI files will dictate the order of the genomes in the ribbon plots.
This script will generate two TSV files: `{prefix}.links.tsv`  `{prefix}.sequence_lengths.tsv`

2. Run the R script
* Required R packages: argparse, gggenomes, gtools, scales
```
usage: plot_synteny_blocks_gggenomes.R [-h] -s SEQUENCES -l LINKS [--scale SCALE] [-p PREFIX]

Plot the ntSynt synteny blocks using gggenomes

optional arguments:
  -h, --help            show this help message and exit
  -s SEQUENCES, --sequences SEQUENCES
                        Input sequence lengths TSV
  -l LINKS, --links LINKS
                        Synteny block links
  --scale SCALE         Length of scale bar in bases (default 1 Gbp)
  -p PREFIX, --prefix PREFIX
                        Output prefix for PNG image (default
                        synteny_gggenomes_plot)
```
Example:
![Example_gggenomes](https://github.com/BirolLab/ntSynt/blob/main/visualization_scripts/example_gggenomes.png)

* These plots are highly customizable, so edit/adapt the script as needed!

## Chromosome sequence painting plots
ggplot2 is used to generate chromosome painting plots to visualize the multi-genome synteny blocks. Here, the synteny blocks for the other species are compared to a selected 'target' species. The segments based on the coordinate system of the 'target' species are coloured based on the chromosome of the other species. In addition, the relative orientation between the 'target' species and the other species is indicated by nudging the boxes up (forward) or down (reverse).

### Steps:
1. Format the ntSynt synteny blocks for the chromosome painting R script using the provided script.
```
usage: format_blocks_chromosome_painting.py [-h] [--convert CONVERT] --target TARGET synteny_tsv

Formatting synteny blocks for chromosome painting

positional arguments:
  synteny_tsv        ntSynt-formatted synteny blocks TSV

optional arguments:
  -h, --help         show this help message and exit
  --convert CONVERT  TSV file with desired conversions for assembly names (optional)
  --target TARGET    Target assembly name
```

This script will output the results in standard out - pipe the output to a file.

2. Run the R script
* Required R packages: argparse, ggplot2, dplyr, tidyr, gtools, scales
```
usage: plot_synteny_blocks-chromosome-painting.R [-h] -b BLOCKS [-g GAPS] [-p PREFIX]

Generate chromosome painting plots to visualize ntSynt synteny

optional arguments:
  -h, --help            show this help message and exit
  -b BLOCKS, --blocks BLOCKS
                        ntSynt-formatted synteny blocks
  -g GAPS, --gaps GAPS  TSV file with gap coordinates. Required headers:
                        chrom, chromStart, chromEnd (optional)
  -p PREFIX, --prefix PREFIX
                        Output prefix for PNG image (optional, default
                        synteny_chromosome)
```
* These plots can also be customized as needed (plot size, scale, comparing additional tools, etc.)

Example:
![Example_gggenomes](https://github.com/BirolLab/ntSynt/blob/main/visualization_scripts/example_chromosome-painting.png)
