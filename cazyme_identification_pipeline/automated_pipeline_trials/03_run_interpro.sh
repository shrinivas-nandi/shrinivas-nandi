#!/usr/bin/env bash
set -euo pipefail

DOMAINS=$1
PROTEOME=$2
RESOLVED=$3

SEQKIT=/home/timothy/programs/seqkit_v2.3.1/seqkit
INTERPRO=/scratch/shrinivas/programs/my_interproscan/interproscan-5.77-108.0/interproscan.sh

mkdir -p output

################################
# Step 1: Clean resolved domain file
# remove .hmm from column 1
################################

awk -F'\t' 'BEGIN{OFS="\t"} {sub(/\.hmm$/, "", $1); print}' \
"$RESOLVED" > output/resolved_domains_clean.tsv

################################
# Step 2: Subset domains of interest
################################

awk -F'\t' '
NR==FNR {d[$1]; next}
{
    for (k in d) {
        if ($1 ~ "^"k"(_|$)") {
            print
            break
        }
    }
}
' "$DOMAINS" output/resolved_domains_clean.tsv \
> output/domains_subset.tsv

################################
# Step 3: Extract protein names
################################

awk -F'\t' '{print $4}' output/domains_subset.tsv | sort -u \
> output/protein_names.txt

################################
# Step 4: Extract sequences
################################

"$SEQKIT" grep -f output/protein_names.txt "$PROTEOME" \
> output/proteins_of_interest.fasta

################################
# Step 5: Run InterProScan
################################

"$INTERPRO" \
  -i output/proteins_of_interest.fasta \
  -f tsv \
  -b output/interpro_annotations \
  -dp \
  -appl Pfam,TIGRFAM, SUPERFAMILY \
  --cpu 60

#,TIGRFAM,SUPERFAMILY,SignalP,TMHMM
# actual run sequence:
# bash extract_domains.sh domains.txt proteome.faa resolved_domains.tsv


