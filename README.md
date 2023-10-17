# regenerationCMap

# Repeating the analysis

## Installation

Clone the repository and use either `devtools::load_all()` or `devtools::install()` to
load the required functions and data. Using `devtools::install()` will also install
package dependencies and is recommended.

## Acquiring required data and preprocessing

"analysis/00.DownloadPreprocess" directory includes code used to compiled required
datasets and some precalculations to speed up the downstream analysis. The output
of the "analysis/00.DownloadPreprocess/00.DownloadPreprocess.R" script is already present
in this repository therefore this script can be skipped. The rest of the scripts should be run
in order to perform the necesarry downloads and pre-processing.


## Running the analysis

Precalculated results for every step can be found within this repository. These are
the instructions to re-generate them from scratch.

'analysis/01.CmapEnrichment' directory contains the code required for running
basic cmap enrichment. Running "01.ScoreCalc.R" will populate the "results" directory
'combined.tsv' in this directory includes enrichment score for comparisons.

'analysis/02.L1000Analysis' contains the code required for running comparisons of results using the old
CMAP dataset with different forms of the new L1000 dataset. Running these scripts in order
will populate respectively named results directories and generate the comparison [heatmap](analysis/02.L1000Analysis/results/cor_plot.png).

