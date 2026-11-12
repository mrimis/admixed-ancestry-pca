#!/usr/bin/env bash
# 00_download.sh
# Download the 1000 Genomes Affymetrix 6.0 array genotypes (GRCh37, 3,450 samples)
# and the sample panel that gives each sample's population label.
# Size: about 750 MB. Run from the repository root.
set -euo pipefail

mkdir -p data
BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip

# -C - lets you resume if the connection drops.
curl -L -C - -o data/affy6.vcf.gz \
  "$BASE/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz"

curl -L -o data/affy6_samples.panel \
  "$BASE/affy_samples.20141118.panel"

echo "Done. Files are in data/"
