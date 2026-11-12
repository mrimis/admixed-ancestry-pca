#!/usr/bin/env bash
# 01_import_qc.sh
# Pick the samples, import the VCF, and run basic sample and variant QC.
# Run from the repository root after 00_download.sh.
set -euo pipefail

mkdir -p work
PANEL=data/affy6_samples.panel

# ---------------------------------------------------------------------------
# 1. Sample lists
#    Cohort    = ACB + ASW  (the admixed samples we pretend are "our study")
#    Reference = YRI (African), CEU (European), CHB (East Asian)
#    CHB is there as a check. It lets us see ancestry that a two-way
#    African/European model would miss.
# ---------------------------------------------------------------------------
pop_ids () { awk -v p="$1" 'NR > 1 && $2 == p { print $1 }' "$PANEL"; }

{ pop_ids ACB; pop_ids ASW; }               > work/cohort.ids
{ pop_ids YRI; pop_ids CEU; pop_ids CHB; }  > work/ref.ids
pop_ids YRI > work/ref_YRI.ids
pop_ids CEU > work/ref_CEU.ids
pop_ids CHB > work/ref_CHB.ids
cat work/cohort.ids work/ref.ids > work/all.ids

# ---------------------------------------------------------------------------
# 2. Import
#    Keep autosomal, biallelic SNPs. Give unnamed variants an ID and drop
#    any variant whose ID appears more than once.
# ---------------------------------------------------------------------------
plink2 --vcf data/affy6.vcf.gz \
  --keep work/all.ids \
  --chr 1-22 \
  --snps-only just-acgt --max-alleles 2 \
  --set-missing-var-ids '@:#:$r:$a' \
  --rm-dup exclude-all \
  --make-pgen --out work/raw

# ---------------------------------------------------------------------------
# 3. QC
#    Variants first, then samples.
#    Note what is NOT here: a Hardy-Weinberg filter on the pooled data.
#    See step 4 for why.
# ---------------------------------------------------------------------------
plink2 --pfile work/raw --geno 0.02 --maf 0.05 --make-pgen --out work/qc_variants
plink2 --pfile work/qc_variants --mind 0.02 --make-pgen --out work/qc

# ---------------------------------------------------------------------------
# 4. Demonstration: why a pooled HWE filter is a mistake here
#    Mixing groups with different allele frequencies creates a deficit of
#    heterozygotes (the Wahlund effect). Good SNPs then fail HWE. The SNPs
#    that fail are the ones that differ most between ancestries, which are
#    the ones PCA needs.
# ---------------------------------------------------------------------------
plink2 --pfile work/qc --hardy --out work/hwe_pooled
plink2 --pfile work/qc --keep work/ref_YRI.ids --hardy --out work/hwe_YRI

count_fail () { awk 'NR > 1 && $NF < 1e-6 { n++ } END { print n + 0 }' "$1"; }
echo "SNPs with HWE p < 1e-6, pooled samples : $(count_fail work/hwe_pooled.hardy)"
echo "SNPs with HWE p < 1e-6, YRI only       : $(count_fail work/hwe_YRI.hardy)"
