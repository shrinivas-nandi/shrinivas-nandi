#!/usr/bin/env bash
set -euo pipefail

################################
# Usage:
# bash 06_run_structure_esm.sh
#
# Looks for:
# /scratch/shrinivas/Sargassum/03_meta_analysis/automated_pipeline_trials/output/domain_split/<DOMAIN>/<DOMAIN>_representatives_ge5.fasta
#
# Runs ESMFold on each file and writes outputs to:
# /scratch/shrinivas/Sargassum/03_meta_analysis/automated_pipeline_trials/output/domain_split/<DOMAIN>/esmfold_output
################################

PARENT_DIR="/scratch/shrinivas/Sargassum/03_meta_analysis/automated_pipeline_trials/output/domain_split"
SIF="/scratch/singularity/ESMFold_v2.0.1_cuda_v11.8.0-rev1.sif"
BIND_DIR="/scratch/shrinivas"
MODEL_DIR="/model"

################################
# Basic checks
################################
if [[ ! -d "$PARENT_DIR" ]]; then
    echo "ERROR: parent directory not found: $PARENT_DIR" >&2
    exit 1
fi

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: Singularity image not found: $SIF" >&2
    exit 1
fi

################################
# Find all representative fasta files
################################
find "$PARENT_DIR" -type f -name "*_representatives_ge5.fasta" | sort | while read -r fasta; do
    DOMAIN_DIR=$(dirname "$fasta")
    DOMAIN=$(basename "$DOMAIN_DIR")
    FASTA_BASENAME=$(basename "$fasta")
    OUTDIR="${DOMAIN_DIR}/esmfold_output"
    LOGFILE="${DOMAIN_DIR}/${DOMAIN}_esmfold.log"

    mkdir -p "$OUTDIR"

    if [[ ! -s "$fasta" ]]; then
        echo "Skipping ${DOMAIN}: fasta is empty -> $fasta"
        continue
    fi

    echo "========================================"
    echo "Domain      : $DOMAIN"
    echo "Input fasta : $fasta"
    echo "Output dir  : $OUTDIR"
    echo "Log file    : $LOGFILE"
    echo "========================================"

    # Optional host-side sanity check
    ls -lh "$fasta"

    singularity exec \
      -B "$BIND_DIR" \
      "$SIF" \
      esm-fold \
      --fasta "$fasta" \
      -o "$OUTDIR" \
      -m "$MODEL_DIR" \
      --cpu-only \
      > "$LOGFILE" 2>&1

    echo "Finished ${DOMAIN}"
done
