#!/usr/bin/env bash
# 03_pca.sh
# Two PCAs, for two different jobs.
#   A. Reference-anchored PCA: describes ancestry. PCs come from the reference
#      samples. The cohort is projected onto them.
#   B. In-cohort PCA: gives covariates for GWAS. PCs come from the unrelated
#      cohort samples. Relatives are projected onto them.
# Run from the repository root after 02_relatedness_ld.sh.
set -euo pipefail

# This follows the projection recipe in the PLINK 2 documentation
# (www.cog-genomics.org/plink/2.0/score, "PCA projection with --score").
# With vcols=chrom,ref,alt the .eigenvec.allele columns are:
#   1 #CHROM  2 ID  3 REF  4 ALT  5 A1  6..15 PC1..PC10
# So variant ID is column 2, the scored allele is column 5, PCs are 6-15.
project () {  # usage: project PCA_PREFIX SAMPLE_LIST OUT_PREFIX
  plink2 --pfile work/qc \
    --keep "$2" \
    --read-freq "$1.acount" \
    --score "$1.eigenvec.allele" 2 5 header-read no-mean-imputation variance-standardize \
    --score-col-nums 6-15 \
    --out "$3"
}

# ---------------------------------------------------------------------------
# A. Reference-anchored PCA
#    Compute on unrelated reference samples only.
#    Then score EVERY sample, reference included, with the same command.
#    That puts reference and cohort on the same scale.
# ---------------------------------------------------------------------------
plink2 --pfile work/qc \
  --keep work/ref_unrel.ids \
  --extract work/prune.prune.in \
  --maf 0.01 \
  --freq counts \
  --pca 10 allele-wts vcols=chrom,ref,alt \
  --out work/ref_pca

project work/ref_pca work/all.ids work/ref_proj
# -> work/ref_proj.sscore

# ---------------------------------------------------------------------------
# B. In-cohort PCA
#    Compute on unrelated cohort samples only.
#    Then score the whole cohort, so relatives get PCs without shaping them.
# ---------------------------------------------------------------------------
plink2 --pfile work/qc \
  --keep work/cohort_unrel.ids \
  --extract work/prune.prune.in \
  --maf 0.01 \
  --freq counts \
  --pca 10 allele-wts vcols=chrom,ref,alt \
  --out work/cohort_pca

project work/cohort_pca work/cohort.ids work/cohort_proj
# -> work/cohort_proj.sscore

# ---------------------------------------------------------------------------
# C. Two PCAs done the wrong way, for comparison in the plots
#    C1: relatives left in.
#    C2: no LD pruning and no removal of long-range LD regions.
# ---------------------------------------------------------------------------
plink2 --pfile work/qc \
  --keep work/cohort.ids \
  --extract work/prune.prune.in \
  --maf 0.01 \
  --pca 10 \
  --out work/wrong_relatives

plink2 --pfile work/qc \
  --keep work/cohort_unrel.ids \
  --maf 0.01 \
  --pca 10 \
  --out work/wrong_noprune

# ---------------------------------------------------------------------------
# D. Inputs for the ancestry estimates in 04_ancestry.R
#    Allele frequencies in each reference group, and cohort genotypes as
#    0/1/2 counts, both on the pruned variant set.
# ---------------------------------------------------------------------------
cut -f2 work/ref_pca.eigenvec.allele | tail -n +2 | sort -u > work/ancestry_snps.txt

for pop in YRI CEU CHB; do
  plink2 --pfile work/qc \
    --keep work/ref_${pop}_unrel.ids \
    --extract work/ancestry_snps.txt \
    --freq --out work/freq_${pop}
done

plink2 --pfile work/qc \
  --keep work/cohort.ids \
  --extract work/ancestry_snps.txt \
  --export A --out work/cohort_geno
# -> work/cohort_geno.raw
