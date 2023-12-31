#!/bin/bash

# Set the source directory for merged files
source_dir="/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication/tree_processing"

# Set the destination directory for aligned files
destination_dir="/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication/tree_processing_aligned"

# Ensure the destination directory exists
mkdir -p "$destination_dir"

# Set the number of threads to use
threads=40  # Change this number based on your system's capabilities

# Run mafft for each input file in parallel
find "$source_dir" -type f -name "merged_tree_ready_*.fa" | xargs -I {} -P "$threads" /home/shrinivas/Programs/mafft-linux64/mafft.bat --localpair --maxiterate 1000 {} ">" "$destination_dir/trial_{}"_aligned.txt


