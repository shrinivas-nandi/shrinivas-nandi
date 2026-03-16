#!/usr/bin/env bash
set -euo pipefail

################################
# Usage:
# bash 05_prepare_domain_dirs.sh domains.txt proteome.faa merged_domain_architecture_filtered.tsv
################################

DOMAINS=$1
PROTEOME=$2
ARCHITECTURE=$3

SEQKIT=/home/timothy/programs/seqkit_v2.3.1/seqkit
OUTDIR="output/domain_split"

mkdir -p "$OUTDIR"

# Activate mmseqs2 once
source ~/mambaforge/etc/profile.d/conda.sh
conda activate mmseqs2

# Master summary
echo -e "domain\tn_proteins\tn_representatives\tn_clusters_ge5" > "${OUTDIR}/cluster_stats.tsv"

################################
# Loop through each domain of interest
################################
while IFS= read -r domain; do
    [[ -z "$domain" ]] && continue

    DOMAIN_DIR="${OUTDIR}/${domain}"
    mkdir -p "$DOMAIN_DIR"

    ################################
    # Extract unique protein headers for this domain
    ################################
    awk -F'\t' -v d="$domain" '
{
    split($3,a,"_")
    if (a[1]==d && $4=="dbCAN") print $1
}' "$ARCHITECTURE"

    ################################
    # Extract matching sequences
    ################################
    "$SEQKIT" grep -f "${DOMAIN_DIR}/${domain}_protein_names.txt" "$PROTEOME" \
        > "${DOMAIN_DIR}/${domain}.fasta"

    n_proteins=$(grep -c '^>' "${DOMAIN_DIR}/${domain}.fasta" || true)
    n_representatives=0
    n_clusters_ge5=0

    ################################
    # Run mmseqs2
    ################################
    if [[ -s "${DOMAIN_DIR}/${domain}.fasta" ]] && [[ "$n_proteins" -gt 0 ]]; then
        (
            cd "$DOMAIN_DIR"

            mmseqs createdb "${domain}.fasta" "${domain}_db"

            mmseqs cluster "${domain}_db" "${domain}_clu" "${domain}_tmp" \
              --min-seq-id 0.30 \
              -c 0.80 \
              --cov-mode 0 \
              --cluster-mode 2

            mmseqs createtsv "${domain}_db" "${domain}_db" "${domain}_clu" "${domain}_clusters.tsv"

            ################################
            # Cluster sizes
            ################################
            awk -F'\t' '
            {
                count[$1]++
            }
            END {
                print "representative\tcluster_size"
                for (rep in count) {
                    print rep "\t" count[rep]
                }
            }' "${domain}_clusters.tsv" | sort -k2,2nr > "${domain}_cluster_sizes.tsv"

            ################################
            # Keep only clusters with >=5 proteins
            ################################
            awk -F'\t' 'NR>1 && $2>=5 {print $1}' "${domain}_cluster_sizes.tsv" \
                > "${domain}_rep_ge5.txt"

            awk -F'\t' 'NR==FNR {keep[$1]=1; next} ($1 in keep)' \
                "${domain}_rep_ge5.txt" "${domain}_clusters.tsv" \
                > "${domain}_clusters_ge5.tsv"

            ################################
            # Representative sequences
            ################################
            mmseqs result2repseq "${domain}_db" "${domain}_clu" "${domain}_repseq_db"
            mmseqs convert2fasta "${domain}_repseq_db" "${domain}_representatives.fasta"

            ################################
            # Extract only representative sequences from clusters >=5
            ################################
            if [[ -s "${domain}_rep_ge5.txt" ]]; then
                "$SEQKIT" grep -f "${domain}_rep_ge5.txt" "${domain}_representatives.fasta" \
                    > "${domain}_representatives_ge5.fasta"
            else
                : > "${domain}_representatives_ge5.fasta"
            fi

            ################################
            # Remove all MMseqs intermediate files
            ################################
            rm -rf \
                "${domain}_db"* \
                "${domain}_clu"* \
                "${domain}_tmp"* \
                "${domain}_repseq_db"*
        )

        n_representatives=$(grep -c '^>' "${DOMAIN_DIR}/${domain}_representatives.fasta" || true)
        n_clusters_ge5=$(grep -c '.' "${DOMAIN_DIR}/${domain}_rep_ge5.txt" || true)
    else
        echo "Skipping ${domain}: no sequences found"
        : > "${DOMAIN_DIR}/${domain}_clusters.tsv"
        : > "${DOMAIN_DIR}/${domain}_cluster_sizes.tsv"
        : > "${DOMAIN_DIR}/${domain}_rep_ge5.txt"
        : > "${DOMAIN_DIR}/${domain}_clusters_ge5.tsv"
        : > "${DOMAIN_DIR}/${domain}_representatives.fasta"
        : > "${DOMAIN_DIR}/${domain}_representatives_ge5.fasta"
    fi

    ################################
    # Per-domain stats
    ################################
    {
        echo -e "domain\tn_proteins\tn_representatives\tn_clusters_ge5"
        echo -e "${domain}\t${n_proteins}\t${n_representatives}\t${n_clusters_ge5}"
    } > "${DOMAIN_DIR}/${domain}_cluster_stats.tsv"

    ################################
    # Append to merged stats
    ################################
    echo -e "${domain}\t${n_proteins}\t${n_representatives}\t${n_clusters_ge5}" >> "${OUTDIR}/cluster_stats.tsv"

done < "$DOMAINS"
