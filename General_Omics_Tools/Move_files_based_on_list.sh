#!/bin/bash

# Source directory containing all files with .fa extension
source_directory="/scratch/shrinivas/Symbiodinium_Meta_Analysis/Cleaned_Reads/Camp_2020"

# Destination directory to move files to
destination_directory="/scratch/shrinivas/Symbiodinium_Meta_Analysis/Cleaned_Reads/Camp_2020/Breviolum"

# Txt file containing the list of filenames to move
files_to_move_txt="SRA_B.txt"

# Ensure destination directory exists
mkdir -p "$destination_directory"

# Move files listed in the txt file to the destination directory
while IFS= read -r partial_string; do
    # Loop over files in source directory
    for file in "$source_directory"/*.fa; do
        # Check if filename contains the partial string
        if [[ "$file" == *"$partial_string"* ]]; then
            # Move the file to the destination directory
            mv "$file" "$destination_directory/"
            echo "Moved: $file"
        fi
    done
done < "$files_to_move_txt"
