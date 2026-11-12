#!/usr/bin/env bash
# 02_relatedness_ld.sh
# Find relatives, build an unrelated sample set, and LD-prune the variants.
# Run from the repository root after 01_import_qc.sh.
set -euo pipefail

# ---------------------------------------------------------------------------
# 1. Relatedness with KING-robust
#    0.0884 is the usual cutoff between second- and third-degree relatives.
#    The table is for you to read. The cutoff run writes the sample lists.
# ---------------------------------------------------------------------------
plink2 --pfile work/qc --make-king-table --king-table-filter 0.0884 --out work/king
plink2 --pfile work/qc --king-cutoff 0.0884 --out work/unrel
# -> work/unrel.king.cutoff.in.id   (keep these)
# -> work/unrel.king.cutoff.out.id  (one member of each related pair)

# Unrelated subsets of the cohort and of each reference group.
keep_unrel () { grep -F -x -f work/unrel.king.cutoff.in.id "$1" > "$2"; }
keep_unrel work/cohort.ids   work/cohort_unrel.ids
keep_unrel work/ref.ids      work/ref_unrel.ids
keep_unrel work/ref_YRI.ids  work/ref_YRI_unrel.ids
keep_unrel work/ref_CEU.ids  work/ref_CEU_unrel.ids
keep_unrel work/ref_CHB.ids  work/ref_CHB_unrel.ids

echo "Cohort:    $(wc -l < work/cohort.ids) samples, $(wc -l < work/cohort_unrel.ids) unrelated"
echo "Reference: $(wc -l < work/ref.ids) samples, $(wc -l < work/ref_unrel.ids) unrelated"

# ---------------------------------------------------------------------------
# 2. LD pruning
#    First drop the known long-range LD regions (GRCh37 coordinates).
#    Then prune what is left. Use unrelated samples to estimate LD.
# ---------------------------------------------------------------------------
plink2 --pfile work/qc \
  --keep work/unrel.king.cutoff.in.id \
  --exclude bed1 resources/high_ld_regions_grch37.txt \
  --indep-pairwise 500kb 0.1 \
  --out work/prune
# -> work/prune.prune.in  (the variant list used from here on)

echo "Variants after QC:      $(grep -vc '^#' work/qc.pvar)"
echo "Variants after pruning: $(wc -l < work/prune.prune.in)"
