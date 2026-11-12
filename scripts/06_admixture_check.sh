#!/usr/bin/env bash
# 06_admixture_check.sh  (optional)
# Cross-check the ancestry estimates from step 4 with ADMIXTURE in supervised mode.
# Needs the `admixture` program on your PATH.
# Run from the repository root after 04_ancestry.R.
set -euo pipefail

command -v admixture > /dev/null || { echo "admixture was not found on your PATH"; exit 1; }

# ---------------------------------------------------------------------------
# 1. Input files
#    ADMIXTURE assumes unrelated samples and reads PLINK 1 .bed files.
#    Use the unrelated reference and cohort samples on the pruned variants.
# ---------------------------------------------------------------------------
plink2 --pfile work/qc \
  --keep work/unrel.king.cutoff.in.id \
  --extract work/prune.prune.in \
  --make-bed --out work/admix

# The .pop file tells ADMIXTURE which samples are references.
# One line per sample, in .fam order. A reference sample gets its group name.
# A cohort sample gets "-", which means "estimate this one".
awk 'NR == FNR { pop[$1] = $2; next }
     { p = pop[$2]; print (p == "YRI" || p == "CEU" || p == "CHB") ? p : "-" }' \
  data/affy6_samples.panel work/admix.fam > work/admix.pop

# ---------------------------------------------------------------------------
# 2. Run ADMIXTURE with K = 3
#    It writes output to the current folder, so run it inside work/.
# ---------------------------------------------------------------------------
( cd work && admixture --supervised admix.bed 3 -j2 > admix.admixture.log )
# -> work/admix.3.Q  ancestry proportions, one row per sample, in .fam order

# ---------------------------------------------------------------------------
# 3. Compare with step 4
# ---------------------------------------------------------------------------
Rscript scripts/06_admixture_check.R
