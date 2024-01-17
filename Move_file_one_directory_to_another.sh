#!/bin/bash

# Source directory containing all files
source_directory="/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication_2/Orthogroup_Sequences"

# Destination directory to move files to
destination_directory="/scratch/shrinivas/Paulinella_Metagenome/Gene_duplication_2/non-target_OG"

# Txt file containing the list of filenames to move
files_to_move_txt="Non_target_orthogroup.txt"

# Ensure destination directory exists
mkdir -p "$destination_directory"

# Move files listed in the txt file to the destination directory
while IFS= read -r filename; do
  if [ -e "$source_directory/$filename" ]; then
    mv "$source_directory/$filename" "$destination_directory/"
    echo "Moved: $filename"
  else
    echo "File not found: $filename"
  fi
done < "$files_to_move_txt"


