#!/usr/bin/env bash
set -euo pipefail

################################
# Usage:
# bash 06_sequence_similarity.sh over100.tsv /path/to/original_hmmer.domtblout
#
# Example:
# bash 06_sequence_similarity.sh over100.tsv /scratch/shrinivas/Sargassum/03_meta_analysis/cazyme_analysis/use_this_domain_resolved_coverage_domtblout.tsv
################################

OVER100=$1
HMMER=$2

BASE_OUTDIR="output/domain_split"
FINAL_OUTDIR="/scratch/shrinivas/Sargassum/03_meta_analysis/cazyme_analysis/output"
SEQKIT=/home/timothy/programs/seqkit_v2.3.1/seqkit
THREADS=60

mkdir -p "$FINAL_OUTDIR"

################################
# Activate mmseqs2
################################
source ~/mambaforge/etc/profile.d/conda.sh
conda activate mmseqs2

################################
# Loop through each selected domain
# Uses first column of over100.tsv
################################
awk 'NF>0 {print $1}' "$OVER100" | while IFS= read -r domain; do
    [[ -z "$domain" ]] && continue

    DOMAIN_DIR="${BASE_OUTDIR}/${domain}"
    SIM_DIR="${DOMAIN_DIR}/sequence_similarity"

    echo "======================================"
    echo "Processing ${domain}"
    echo "======================================"

    if [[ ! -d "$DOMAIN_DIR" ]]; then
        echo "Skipping ${domain}: ${DOMAIN_DIR} not found"
        continue
    fi

    REP_TXT="${DOMAIN_DIR}/${domain}_rep_ge5.txt"
    REP_FASTA="${DOMAIN_DIR}/${domain}_representatives_ge5.fasta"
    CLUSTERS_GE5="${DOMAIN_DIR}/${domain}_clusters_ge5.tsv"

    if [[ ! -s "$REP_TXT" ]]; then
        echo "Skipping ${domain}: missing or empty ${REP_TXT}"
        continue
    fi

    if [[ ! -s "$REP_FASTA" ]]; then
        echo "Skipping ${domain}: missing or empty ${REP_FASTA}"
        continue
    fi

    ################################
    # 1. Clean old noise from domain dir
    ################################
    rm -f "${DOMAIN_DIR}/${domain}_clu."*
    rm -f "${DOMAIN_DIR}/${domain}_db"*
    rm -rf "${DOMAIN_DIR}/${domain}_tmp"*

    ################################
    # 2. Create fresh sequence similarity dir
    ################################
    mkdir -p "$SIM_DIR"

    ################################
    # 3. Copy ge5 details into similarity dir
    ################################
    cp "$REP_TXT" "${SIM_DIR}/"
    cp "$REP_FASTA" "${SIM_DIR}/"
    [[ -f "$CLUSTERS_GE5" ]] && cp "$CLUSTERS_GE5" "${SIM_DIR}/"

    (
        cd "$SIM_DIR"

        REP_TXT_LOCAL="${domain}_rep_ge5.txt"
        REP_FASTA_LOCAL="${domain}_representatives_ge5.fasta"

        DOMAIN_COORDS_TSV="${domain}_domain_region_coordinates.tsv"
        DOMAIN_COORDS_BED="${domain}_domain_region_coordinates.bed"
        DOMAIN_ONLY_FASTA="${domain}_domain_only_reps_ge5.fasta"

        FULL_DB="${domain}_representatives_ge5_DB"
        FULL_RES="${domain}_representatives_ge5_allvsall"
        FULL_TMP="tmp_${domain}_representatives_ge5"
        FULL_TSV="${domain}_representatives_ge5_allvsall.tsv"

        DOM_DB="${domain}_domain_only_reps_ge5_DB"
        DOM_RES="${domain}_domain_only_reps_ge5_allvsall"
        DOM_TMP="tmp_${domain}_domain_only_reps_ge5"
        DOM_TSV="${domain}_domain_only_reps_ge5_allvsall.tsv"

        ################################
        # 4. Extract domain coordinates from original HMMER file
        ################################
        awk -v d="$domain" 'BEGIN{FS=OFS="\t"}
        NR==FNR {keep[$1]=1; next}
        $1==(d ".hmm") && ($4 in keep) {
            print $4, $18, $19
        }' "$REP_TXT_LOCAL" "$HMMER" > "$DOMAIN_COORDS_TSV"

        if [[ -s "$DOMAIN_COORDS_TSV" ]]; then
            ################################
            # 5. Convert to BED
            ################################
            awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3}' "$DOMAIN_COORDS_TSV" > "$DOMAIN_COORDS_BED"

            ################################
            # 6. Extract domain-only sequences
            ################################
            "$SEQKIT" subseq \
              --bed "$DOMAIN_COORDS_BED" \
              "$REP_FASTA_LOCAL" \
              > "$DOMAIN_ONLY_FASTA"
        else
            echo "No coordinates found for ${domain}"
            : > "$DOMAIN_ONLY_FASTA"
        fi

        ################################
        # 7. Full-length all-vs-all
        ################################
        if [[ -s "$REP_FASTA_LOCAL" ]]; then
            mmseqs createdb "$REP_FASTA_LOCAL" "$FULL_DB"

            mmseqs search \
              "$FULL_DB" \
              "$FULL_DB" \
              "$FULL_RES" \
              "$FULL_TMP" \
              --threads "$THREADS" \
              -s 7.5 \
              -e 1e-3 \
              --min-seq-id 0.0 \
              --cov-mode 0

            mmseqs convertalis \
              "$FULL_DB" \
              "$FULL_DB" \
              "$FULL_RES" \
              "$FULL_TSV" \
              --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov"
        else
            : > "$FULL_TSV"
        fi

        ################################
        # 8. Domain-only all-vs-all
        ################################
        if [[ -s "$DOMAIN_ONLY_FASTA" ]]; then
            mmseqs createdb "$DOMAIN_ONLY_FASTA" "$DOM_DB"

            mmseqs search \
              "$DOM_DB" \
              "$DOM_DB" \
              "$DOM_RES" \
              "$DOM_TMP" \
              --threads "$THREADS" \
              -s 7.5 \
              -e 1e-3 \
              --min-seq-id 0.0 \
              --cov-mode 0

            mmseqs convertalis \
              "$DOM_DB" \
              "$DOM_DB" \
              "$DOM_RES" \
              "$DOM_TSV" \
              --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov"
        else
            : > "$DOM_TSV"
        fi

        ################################
        # 9. Clean MMseqs/temp noise
        ################################
        rm -rf "$FULL_TMP" "$DOM_TMP"

        rm -f "${FULL_RES}."[0-9]* "${FULL_RES}.dbtype" "${FULL_RES}.index"
        rm -f "${DOM_RES}."[0-9]* "${DOM_RES}.dbtype" "${DOM_RES}.index"

        rm -f "${FULL_DB}"*
        rm -f "${DOM_DB}"*

        rm -f "$DOMAIN_COORDS_TSV" "$DOMAIN_COORDS_BED"
        rm -f "${REP_FASTA_LOCAL}.seqkit.fai"
    )

    ################################
    # 10. Copy final similarity outputs back to main output dir
    ################################
    if [[ -s "${SIM_DIR}/${domain}_representatives_ge5_allvsall.tsv" ]]; then
        cp "${SIM_DIR}/${domain}_representatives_ge5_allvsall.tsv" \
           "${FINAL_OUTDIR}/${domain}_representatives_ge5_allvsall.tsv"
    fi

    if [[ -s "${SIM_DIR}/${domain}_domain_only_reps_ge5_allvsall.tsv" ]]; then
        cp "${SIM_DIR}/${domain}_domain_only_reps_ge5_allvsall.tsv" \
           "${FINAL_OUTDIR}/${domain}_domain_only_reps_ge5_allvsall.tsv"
    fi

done
