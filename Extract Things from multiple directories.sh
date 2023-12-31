Extract Things from multiple directories.sh

find /scratch/shrinivas/Paulinella_Metagenome/PaulHet_DNA.spades.scaffolds.fasta.MetaBAT2_bins/busco_output_correct -type f -name "*.txt" -exec mv {} /scratch/shrinivas/Paulinella_Metagenome/PaulHet_DNA.spades.scaffolds.fasta.MetaBAT2_bins/busco_output_correct/summary_all \;


#!/bin/bash

# Change to the summary_all directory
cd /scratch/shrinivas/Paulinella_Metagenome/PaulHet_DNA.spades.scaffolds.fasta.MetaBAT2_bins/busco_output_correct/summary_all

# Create the "all_results.txt" file if it doesn't exist
touch all_results.txt

# Iterate through all .txt files in the directory and its subdirectories
find . -type f -name "*.txt" | while read -r file; do
    # Get the file name without the path
    filename=$(basename "$file")

    # Add a subheading with the filename to "all_results.txt"
    echo "=== $filename ===" >> all_results.txt

    # Use awk to print lines after "***** Results: *****" and append them to "all_results.txt"
    awk '/\*\*\*\*\* Results: \*\*\*\*\*/ {p=1; next} p' "$file" >> all_results.txt
done
