#!/usr/bin/env bash
# run_all.sh : run the whole tutorial from a clean start.
# Needs plink2 and R (with data.table and ggplot2) on your PATH.
set -euo pipefail
[ -f data/affy6.vcf.gz ] || bash scripts/00_download.sh
bash scripts/01_import_qc.sh
bash scripts/02_relatedness_ld.sh
bash scripts/03_pca.sh
Rscript scripts/04_ancestry.R
Rscript scripts/05_plots.R
