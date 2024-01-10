In this script we focus on sorting out the downsampling bits of code that are relatively more general. 

Downsample hits after balst result is obtained. In this case we use a method already developed by Tim stephens and Dana pierce. 
This script basically removes hits that belong to the same species so that each spp is represented only once. 

zcat extracted_sequences_blast.tsv | /home/timothy/GitHub/Utils/Blast/taxonomically_downsample_hits/taxonomically_downsample_hits.sh --taxid_col 16 > downsampled_hits.tsv
wk 'NR>1 {split($1, a, "\\."); gsub(/^[ \t]+|[ \t]+$/, "", a[1]); print a[1]}' new_output.tsv | sort -u | while read -r value; do
    grep -E "^${value}" output.tsv > "fasta_headers/${value}.txt"
done

#########


# Create the 'fasta_headers' directory if it doesn't exist
mkdir -p fasta_headers

# Assuming 'all_complete_seq.fasta' is in the same directory
all_seq_file="downsampled_seqs.fa"

# Iterate through each header file
for header_file in fasta_headers/*.txt; do
    # Extract the header value from the file name (remove the path and extension)
    header_value=$(basename "$header_file" .txt)

    # Use grep to find the matching sequence in 'all_complete_seq.fa' and append to a new file
    grep -A 1 -F ">$header_value" "$all_seq_file" >> "$header_file"
done


# lets try it for one file
#!/bin/bash

# Specify file paths
HEADER_FILE="OG0000006.txt"
FASTA_FILE="downsampled_seqs.fa"
OUTPUT_FILE="OG0000006.fasta"

# Use seqkit to extract sequences based on headers
/home/timothy/programs/seqkit_v2.3.1/seqkit grep -f "$HEADER_FILE" "$FASTA_FILE" > "$OUTPUT_FILE"

echo "Extraction completed."
