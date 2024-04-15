e###Workflow for bootstrap tree generation 

##Step 1: Select genomes for each of the species of interest 
#Normally use python or something to do this selection

##Step 2: Utilize orthofinder to generate orthogroups 
/home/timothy/programs/OrthoFinder_v2.5.4/orthofinder -f /scratch/shrinivas/Paulinella_Metagenome/Tree_Analysis/Cyanobacteria_Genomes -t 50

## Step 3: Select orthogroups of interest based on the table generated
#Apply a specific criteria


## Step 4: Move the selected orthogroups to another directory 

source_directory=<dir>
destination_directory=<dir>

while read -r filename; do
    cp "$source_directory/$filename.fa" "$destination_directory/"
done < Orthogroup_List_for_Tree_27_10.txt

## Step 4: Mafft align each orthogroup sequence file 

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


### Step 5: Trimal trimal aligned files 
# Trimal loop
#!/bin/bash

# Set the input directory path
input_dir="/input/directory"

# Set the output directory path
output_dir="enter/output/directory "

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all files in the input directory with a .fasta extension
for file in "$input_dir"/*.fasta; do
    # Extract the filename without the extension
    filename=$(basename "$file" .fasta)

    # Define the output file path in the new directory
    output_file="$output_dir/${filename}_trim.fasta"

    # Run trimal command
    trimal -automated1 -in "$file" -out "$output_file"

    echo "Processed: $filename"
done



## Step 6: Run Iqtree for each individual aligned and trimmed file 
#!/bin/bash
source_directory="/scratch/shrinivas/Paulinella_Metagenome/Oct25_2023/Oct_27/trimal_trimmed"

# Destination directory for IQ-TREE output
destination_directory="/scratch/shrinivas/Paulinella_Metagenome/Oct25_2023/Oct_27/iqtree_results"

# Ensure the destination directory exists
mkdir -p "$destination_directory"

# Loop through each file in the source directory
for source_file in "$source_directory"/*.fa; do
    if [ -f "$source_file" ]; then
        filename=$(basename "$source_file")
        output_filename="${filename%.*}.tree"
        output_filepath="$destination_directory/$output_filename"
        
        # Run the iqtree command with -nt AUTO
        echo "/home/timothy/programs/iqtree-1.6.12-Linux/bin/iqtree -s $source_file -m TEST -bb 1000 -nt AUTO -pre $output_filepath"
        
        echo "Processed: $filename"
    fi
done | parallel -j 10
echo "IQ-TREE process complete."



## Step 7: Extract best fit model: Export results from Iqtree log file for 
#!/bin/bash

# Specify the directory where your .log files are located
log_directory="/scratch/shrinivas/Paulinella_Metagenome/Oct25_2023/Oct_27/iqtree_results"

# Specify the output TSV file
output_file="/scratch/shrinivas/Paulinella_Metagenome/Oct25_2023/Oct_27/best_fit_models.tsv"

# Initialize the output file with a header
echo -e "File Name\tBest-fit Model Line" > "$output_file"

# Loop through all .log files in the directory
for log_file in "$log_directory"/*.log; do
    if [ -f "$log_file" ]; then
        best_fit_model_line=$(grep "Best-fit model" "$log_file")
        if [ -n "$best_fit_model_line" ]; then
            # Extract just the file name from the full path
            file_name=$(basename "$log_file")
            # Append the result to the output file in TSV format
            echo -e "$file_name\t$best_fit_model_line" >> "$output_file"
        fi
    fi
done

Clean up best fit model file 
# Step 7b. Look to check if all of the columns have values and remove any that fail 
#!/bin/bash

# Specify the path to your TSV file
tsv_file="your_file.tsv"

# Flag to indicate if any row is missing values
missing_values=0

# Iterate over each line in the file
while IFS=$'\t' read -r col1 col2 _; do
    # Check if either column is empty
    if [[ -z "$col1" || -z "$col2" ]]; then
        echo "Found row with missing values: $col1 $col2"
        missing_values=1
    fi
done < "$tsv_file"

# Check if any row has missing values
if [ "$missing_values" -eq 0 ]; then
    echo "All rows have values in both columns."
else
    echo "Some rows have missing values."
fi

# remove
#!/bin/bash

# Specify the path to your input and output TSV files
input_tsv="input_file.tsv"
output_tsv="output_file.tsv"

# Remove any existing output file
rm -f "$output_tsv"

# Iterate over each line in the input file
while IFS=$'\t' read -r col1 col2 _; do
    # Check if both columns have values
    if [[ -n "$col1" && -n "$col2" ]]; then
        # Write the row to the output file
        echo -e "$col1\t$col2" >> "$output_tsv"
    fi
done < "$input_tsv"



# Step 8: Generate concatenated.out file and partitions.txt file
## To do this properly, the names of the sequences need to be correct. ie need to be the same so they merge correctly. 
## I do it in steps, by removing bits and pieces in seq headers that are common. For example in this case the sequence headers were like prochloroccoccus___1 and prochlorococcus___2 
## So i delete everything after ___

#!/bin/bash
# Define the directory containing the FASTA files
sequence_directory="/scratch/shrinivas/Paulinella/Paulinella_Metagenome/Tree_Analysis/SAR_tree_20_2_24/tree_building/fasta_for_build"

# Iterate through the FASTA files in the specified directory
for sequence_file in "$sequence_directory"/*.fasta; do
    # Create a new file to store modified sequences
    new_file="${sequence_file%.fasta}_modified.fasta"
    
    # Process each line in the FASTA file
    while IFS= read -r line; do
        # Check if the line is a header line (starts with ">")
        if [[ $line =~ ^\> ]]; then
            # Check if the header contains the string "___"
            if [[ $line == *___* ]]; then
                # Truncate the header at the position of "___"
                truncated_header="${line%%___*}"
                echo "$truncated_header" >> "$new_file"
            # Check if the header starts with ">g"
            elif [[ $line == ">g"* ]]; then
                # Keep only the first two characters of the header
                short_header="${line:0:2}"
                echo "$short_header" >> "$new_file"
            else:
                # Otherwise, truncate the header to the first 8 characters
                short_header=$(echo "$line" | cut -c 1-15)
                echo "$short_header" >> "$new_file"
            fi
        else:
            # If it's not a header line, write it as-is
            echo "$line" >> "$new_file"
        fi
    done < "$sequence_file"
done


AMAS

# Make the partitions and concat out file 
python AMAS.py concat -f fasta -d dna -i *fas -u nexus
python /home/shrinivas/Programs/AMAS/amas/AMAS.py concat -f fasta -d aa -i /scratch/shrinivas/Paulinella_Metagenome/Tree_Analysis/Cyanobacteria_Genomes/aligned_trimmed/* --concat-part trial.partitions.txt --part-format nexus --concat-out trial.concatenated.out --out-format fasta --cores 40



# Step 9: Generate nexu file for bootstrap tree 

# cleaning up file 
sed -i 's/msa.*corrected/msa corrected/g' your_file.txt

## Add best fit model to partition file

'''
# nexus file format
charset OG0116139 = 626031-626207;
charset OG0116585 = 626208-626371;
charset OG0116681 = 626372-627105;
charpartition mine = DCMut+F+I+G4:OG0018548,LG+G4:OG0020202,LG+F+G4:OG0021206,LG+G4:OG0021527,LG+I+G4:OG0021566,LG+F+G4:OG0022230,LG+F+G4:OG0024750,LG+G4:OG0024870,LG+I+G4:O
G0026270,LG+G4:OG0026325,LG+G4:OG0026836,WAG+F+G4:OG0028086,LG+G4:OG0028225,LG+I+G4:OG0028342,LG+I+G4:OG0028380,LG+G4:OG0030367,LG+G4:OG0030431,PMB+I+G4:OG0030486,WAG+G4:OG0
030511,LG+F+I+G4:OG0030942,rtREV+F

'''

# Step 10:Run Bootstrap tree 
iqtree -s concatenated.out -spp partitions.txt -bb 2000 -nt AUTO -ntmax 48


####################TROUBLE SHOOTING######################

# count the number of log files

find /scratch/shrinivas/Paulinella/Paulinella_Metagenome/Tree_Analysis/Cyanobacteria_Genomes/iqtree_results/log_files -type f -name "*.log" | wc -l

# 2236 log files 

# how many models in best fit models
wc -l file.tsv

# 2094 rows 

#number of mafft align files 

find -type f -name "*.fasta" | wc -l
# 2234

# cleaned up the list to remove any missing values 

# make the first one a list

awk -F'\t' '{print $1}' output_file.tsv > candidate_list.txt

## now bring in the files from the other directory 

/scratch/shrinivas/Paulinella/Paulinella_Metagenome/Tree_Analysis/Cyanobacteria_Genomes/concated_tree/candidate_list.txt

# location of fasta files 


#!/bin/bash

# Source directory
source_dir="/scratch/shrinivas/Paulinella/Paulinella_Metagenome/Tree_Analysis/Cyanobacteria_Genomes/aligned_trimmed/phase1/phase2/corrected_headers"

# Destination directory
destination_dir="/scratch/shrinivas/Paulinella/Paulinella_Metagenome/Tree_Analysis/Cyanobacteria_Genomes/concated_tree/fasta_files"

# Iterate over each entry in candidate_list.txt
while IFS= read -r entry; do
    # Extract the first 9 characters from the entry
    candidate="${entry:0:9}"

    # Search for files in the source directory matching the candidate
    matching_files=("$source_dir/$candidate"*)
    
    # Check if any files were found
    if [ ${#matching_files[@]} -gt 0 ]; then
        # Copy matching files to the destination directory
        cp "${matching_files[@]}" "$destination_dir"
    fi
done < candidate_list.


# lets take the new counts
# best_fit_model = 2091
# number of fasta_files = 2091

# now lets clean this shit up 

##### Now go ahead in the actual tsv file and fix it 
sed 's/_msa_align_trim\.fasta\.fasta_iqtree\.log//g' repeat_removed.tsv > new_file.tsv


awk '{gsub(/ *: */, ":", $0); gsub(/ +/, " ", $0); print $0}' 03_tmp.txt > 04_tmp.txt

awk '{print $0 ","}' 04_tmp.txt > 05_tmp.txt

awk '{gsub(/ *: */, ":", $0); gsub(/ +/, " ", $0); print $0}' 05_tmp.txt > 06_tmp.txt

sed 's/ ,/,/g' 06_tmp.txt > 07_tmp.txt

tr '\n' ' ' < 07_tmp.txt > best_fit_models.txt


sed 's/_msa_align_trim_header_fixed_phase1_header_fixed_phase2_corrected//g' partitions.txt > tmp_partitions.txt



sed 's/p.*_/_/g' tmp_partitions.txt > 2_tmp_partitions.txt


sed 's/_//g' 2_tmp_partitions.txt > 3_tmp_partitions.txt


tr '\n' ' ' < best_fit_model.txt | paste -sd ' ' - > merged_best_fit.txt


iqtree -s concatenated.out -spp nexus_file.txt -bb 2000 -nt AUTO -ntmax 48



