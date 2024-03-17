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




# Step 8: Generate concatenated.out file and partitions.txt file
## To do this properly, the names of the sequences need to be correct. ie need to be the same so they merge correctly. 
## I do it in steps, by removing bits and pieces in seq headers that are common. For example in this case the sequence headers were like prochloroccoccus___1 and prochlorococcus___2 
## So i delete everything after ___

#!/bin/bash

# Directory containing the sequence files
sequence_directory="/scratch/shrinivas/Paulinella_Metagenome/Tree_Analysis/Cyanobacteria_Genomes/aligned_trimmed"

# Iterate through sequence files
for sequence_file in "$sequence_directory"/*.fasta; do
    # Create a new file to store modified sequences with the "_header_fixed_phase1" suffix
    new_file="${sequence_file%.fasta}_header_fixed_phase1.fasta"
    
    # Write modified sequences to the new file
    while IFS= read -r line; do
        # Check if the line is a sequence header
        if [[ $line =~ ^\> ]]; then
            # Find the position of the third "_" in the header
            position=$(echo "$line" | grep -b -o '___' | cut -d: -f1 | head -n 1)
            
            # If "___" is found, truncate the header
            if [ "$position" ]; then
                truncated_header="${line:0:$position-1}"
                echo "$truncated_header" >> "$new_file"
            else
                echo "$line" >> "$new_file"
            fi
        else
            # If it's not a header, simply write the line to the new file
            echo "$line" >> "$new_file"
        fi
    done < "$sequence_file"
done

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


