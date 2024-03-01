## This script is desgined for blasting, taking top hits, downsampling, retrieving sequences and then building a phylogenetic tree for all candidates. 
## If building a regular phylogenetic tree just follow the last three steps 


# Step 0: Orthofinder analysis ################################################################################################################################
/home/timothy/programs/OrthoFinder_v2.5.4/orthofinder -f <dir>

# Step 1: Isolate candidate orthogroups after running orthofinder ################################################################
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



# Step 2: Isolate representative sequences ################################################################################################
#!/bin/bash

output_file="extracted_sequences.fa"

# Remove the existing output file if it exists
[ -e "$output_file" ] && rm "$output_file"

for file in /scratch/shrinivas/Paulinella_Metagenome/HGT/candidate_OG/*.fa; do
  /home/timothy/programs/seqkit_v2.3.1/seqkit sort -l -r -w 0 "$file" | head -n 2 | awk -v filename="$(basename "$file")" 'NR==1 {print ">" filename "_" $1} NR==2 {print $0}' >> "$output_file"
done


# Step 3: Perform a blastp ################################################################################################

/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond blastp --ultra-sensitive --iterate --max-target-seqs 0 --evalue 0.00001 --query /scratch/shrinivas/Paulinella_Metagenome/HGT --db /scratch/databases/ncbi/2022_07/nr.dmnd --out Orthogroup_representative_BLAST.tsv --threads 70 --compress 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms


# Step 4: Downsample the blast results 
zcat query.fa.domaind_blastp.outfmt6.gz | /home/timothy/GitHub/Utils/Blast/taxonomically_downsample_hits/taxonomically_downsample_hits.sh --taxid_col 16 --target_taxa class --num_seqs_per_taxid 1 > downsampled_hits.tsv

# Step 5: Extract blast sequences
/scratch/databases/bin/DIAMOND_v2.0.15/bin/diamond getseq --db /scratch/databases/ncbi/2022_07/nr.dmnd | /home/timothy/GitHub/Utils/Blast/taxonomically_downsample_hits/grepf_fasta.py -f <(cut -f2 downsampled_hits.tsv) > downsampled_seqs.fa

# Step 6: Split based on sequence headers (Needs to be tested)
## Create mergeable output files with the OG orthogroup
#!/bin/bash

# Set the directories
dir1="/path/to/dir1"
dir2="/path/to/dir2"
output_dir="/path/to/output_directory"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Iterate over each file in dir1
for file1 in "$dir1"/*; do
    # Get the filename without the directory path
    filename=$(basename "$file1")
    
    # Extract the shared part of the filename
    shared_part="${filename%%_*}"
    
    # Search for corresponding file in dir2
    file2=$(find "$dir2" -type f -name "${shared_part}_*")
    
    if [ -n "$file2" ]; then
        # Merge the contents of the two files into a new file in the output directory
        cat "$file1" "$file2" > "$output_dir/$filename"
        echo "Merged $filename"
    else
        echo "No corresponding file found for $filename in $dir2"
    fi
done

## Step 7: Mafft 
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


## Step 8: Trimal
#!/bin/bash

# Set the directory containing the files
input_directory="/path/to/directory"

# Iterate over each file in the directory
for file in "$input_directory"/*.fasta; do
    # Extract the filename without extension
    filename=$(basename "$file" .fasta)
    
    # Run trimal command
    trimal -automated1 -in "$file" -out "${input_directory}/${filename}_trim.fasta"
done



## Step 9: iqtree



## Step 10: physorter







