#!/bin/bash

dir1="/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication/Phylogenetic_Analysis_Orthogroups"
dir2="/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication/tree_processing/incorrectly_done_files"

# Loop through files in directory1
for file1 in "$dir1"/*; do
    # Extract the last 12 characters of the file1 name
    suffix=$(echo "$file1" | rev | cut -c 1-12 | rev)

    # Search for matching files in directory2
    matching_file2=$(find "$dir2" -type f -name "*$suffix")

    # If a match is found, concatenate the files
    if [ -n "$matching_file2" ]; then
        concat_file="concat_$(basename "$file1")"
        cat "$file1" "$matching_file2" > "$concat_file"
        echo "Concatenated: $concat_file"
    else
        echo "No matching file found for $file1"
    fi
done

:/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication/tree_processing
/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication/down_sampledhits/files_not_yet_run


cat Orthogroup_to_extract.txt | xargs -I {} mv {}.fa /scratch/shrinivas/Paulinella_Metagenome/Oct24_2023/Specific_Paulinella_Only/OrthoFinder/Results_Oct24/Orthogroup_Sequences/upset_analysis

#!/bin/bash
egg
