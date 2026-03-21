#!/usr/bin/env bash
set -euo pipefail

################################
# Usage:
# bash 07_run_structure_esm.sh
#
# Looks for files like:
# /scratch/shrinivas/Sargassum/03_meta_analysis/automated_pipeline_trials/output/domain_split/<DOMAIN>/sequence_similarity/*_representatives_ge5.fasta
# /scratch/shrinivas/Sargassum/03_meta_analysis/automated_pipeline_trials/output/domain_split/<DOMAIN>/sequence_similarity/*_domain_only_reps_ge5.fasta
#
# Runs ESMFold on each file and writes outputs into the same sequence_similarity folder:
# - esmfold_output_representatives/
# - esmfold_output_domain_only/
################################

PARENT_DIR="/scratch/shrinivas/Sargassum/03_meta_analysis/cazyme_analysis/output/domain_split"
SIF="/scratch/singularity/ESMFold_v2.0.1_cuda_v11.8.0-rev1.sif"
BIND_DIR="/scratch/shrinivas"
MODEL_DIR="/model"

################################
# Activate environment
################################
#source ~/mambaforge/etc/profile.d/conda.sh
#conda activate singularity

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
# Find all relevant fasta files in sequence_similarity subdirs
################################
find "$PARENT_DIR" -type f \( -name "*_representatives_ge5.fasta" -o -name "*_domain_only_reps_ge5.fasta" \) | sort | while read -r fasta; do

    SEQSIM_DIR=$(dirname "$fasta")
    PARENT_DOMAIN_DIR=$(dirname "$SEQSIM_DIR")
    DOMAIN=$(basename "$PARENT_DOMAIN_DIR")
    FASTA_BASENAME=$(basename "$fasta")

    # Only process files inside a sequence_similarity directory
    if [[ "$(basename "$SEQSIM_DIR")" != "sequence_similarity" ]]; then
        continue
    fi

    ################################
    # Decide output directory + log name
    ################################
    if [[ "$FASTA_BASENAME" == *_domain_only_reps_ge5.fasta ]]; then
        OUTDIR="${SEQSIM_DIR}/esmfold_output_domain_only"
        PREFIX="${FASTA_BASENAME%_domain_only_reps_ge5.fasta}"
        LOGFILE="${SEQSIM_DIR}/${PREFIX}_domain_only_esmfold.log"
        RUN_TYPE="domain_only"
    elif [[ "$FASTA_BASENAME" == *_representatives_ge5.fasta ]]; then
        OUTDIR="${SEQSIM_DIR}/esmfold_output_representatives"
        PREFIX="${FASTA_BASENAME%_representatives_ge5.fasta}"
        LOGFILE="${SEQSIM_DIR}/${PREFIX}_representatives_esmfold.log"
        RUN_TYPE="representatives"
    else
        echo "Skipping unrecognized file pattern: $fasta"
        continue
    fi

    mkdir -p "$OUTDIR"

    if [[ ! -s "$fasta" ]]; then
        echo "Skipping ${DOMAIN} (${RUN_TYPE}): fasta is empty -> $fasta"
        continue
    fi

    echo "========================================"
    echo "Domain      : $DOMAIN"
    echo "Run type    : $RUN_TYPE"
    echo "Input fasta : $fasta"
    echo "Output dir  : $OUTDIR"
    echo "Log file    : $LOGFILE"
    echo "========================================"

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

    echo "Finished ${DOMAIN} (${RUN_TYPE})"
done
