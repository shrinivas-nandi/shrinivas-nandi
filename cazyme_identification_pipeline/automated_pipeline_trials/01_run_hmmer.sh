#!/usr/bin/env bash
# 01_run_hmmer.sh
set -euo pipefail

################################
# Usage
################################
# bash 01_run_hmmer.sh proteome.fasta

PROTEOME=$1

THREADS=64
COVERAGE=0.60
DBCAN="/scratch/shrinivas/Databases/CaZy_db/dbCAN-HMMdb-V14.txt"

mkdir -p output
source ~/mambaforge/etc/profile.d/conda.sh

################################
# Step 1: HMMER
################################

conda activate hmmer

hmmscan \
  --cpu "$THREADS" \
  --domtblout output/dbcan.domtblout \
  --tblout output/dbcan.tblout \
  -E 1e-15 --domE 1e-15 \
  "$DBCAN" \
  "$PROTEOME" \
  > output/hmmer.log

################################
# Step 2: Convert domtblout to true TSV and append coverage
# output stays headerless
################################

awk 'BEGIN{OFS="\t"}
!/^#/ {
    coverage = ($17 - $16 + 1) / $3
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,coverage
}' output/dbcan.domtblout > output/dbcan_with_cov.tsv

################################
# Step 3: Coverage filter
################################

awk -F'\t' -v c="$COVERAGE" 'BEGIN{OFS="\t"} $23 >= c' \
  output/dbcan_with_cov.tsv > output/dbcan_cov_filtered.tsv
