#!/bin/bash

# Source directory containing all files with .fa extension
source_directory="/scratch/shrinivas/Paulinella_Metagenome/Oct24_2023/Outgroup_included_analysis/OrthoFinder/Results_Oct24_1/Orthogroup_Sequences"

# Destination directory to move files to
destination_directory="/scratch/shrinivas/Paulinella_Metagenome/HGT/candidate_OG"

# Txt file containing the list of filenames to move
files_to_move_txt="files2move.txt"

# Ensure destination directory exists
mkdir -p "$destination_directory"

# Move files listed in the txt file to the destination directory
while IFS= read -r filename; do
    filename_with_extension="$filename.fa"
    if [ -e "$source_directory/$filename_with_extension" ]; then
        mv "$source_directory/$filename_with_extension" "$destination_directory/"
        echo "Moved: $filename_with_extension"
    else
        echo "File not found: $filename_with_extension"
    fi
done < "$files_to_move_txt"

