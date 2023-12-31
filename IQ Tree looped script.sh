#IQ Tree looped script.sh

#!/bin/bash
source_directory="/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication/tree_processing_aligned/trimal_aligned"

# Destination directory for IQ-TREE output
destination_directory="/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication/tree_processing_aligned/trimal_aligned/iqtree_results"

# Ensure the destination directory exists
mkdir -p "$destination_directory"

# Loop through each file in the source directory
for source_file in "$source_directory"/*.fasta; do
    if [ -f "$source_file" ]; then
        filename=$(basename "$source_file")
        output_filename="${filename%.*}.tree"
        output_filepath="$destination_directory/$output_filename"
        
        # Run the iqtree command with -nt AUTO
        /home/timothy/programs/iqtree-1.6.12-Linux/bin/iqtree -s "$source_file" -m TEST -bb 1000 -nt 26 -pre "$output_filepath" 
        
        echo "Processed: $filename"
    fi
done

echo "IQ-TREE process complete."