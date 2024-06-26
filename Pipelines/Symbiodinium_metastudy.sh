Symbiodinium_metastudy.sh

# PHASE 1: GENERATING COUNTS TABLES (ON SERVER)

'''
Downloading the raw reads 
Utilize kingfisher wwood to do this.
For trial puprose we will use this PRJNA342240
'''
# lets get the SRA files associated with this 
module run Kingfisher_v0.4.1-rev1 kingfisher annotate --bioprojects PRJNA342240 -o PRJNA342240.tsv -f tsv
cut -f1 your_file.tsv > PRJ.txt # to get only the run file 

# script to download all SRA associated
#!/bin/bash

# Assuming your list is stored in a file named "input_list.txt"
input_list="input_list.txt"

# Loop through each item in the list
while IFS= read -r line; do
    # Run the module command with the specified arguments
    module run Kingfisher_v0.4.1-rev1 kingfisher get -r "$line" -m ena-ascp aws-http prefetch
done < "$input_list"

# Check qualtiy 
fastqc *.fastq
mutliqc .

# Clean up the fasta files (for single end reads)
fastp -i SRR4237276.fastq -o SRR4237276_clean.fastq --failed_out SRR4237276_failed_out.txt --qualified_quality_phred 20 --unqualified_percent_limit 10 --thread 16 

# (for paired end reads)
fastp -i {file}_1.fastq -I {file}_2.fastq -o {file}_cleaned_1.fastq -O {file}_cleaned_2.fastq  --failed_out SRR4237276_failed_out.txt
--qualified_quality_phred 20 --unqualified_percent_limit 10 --thread 16 --detect_adapter_for_pe -f 20 -t 20

# Now once clean. We run SALMON (quasi mapping. doesnt do gene prediction. Good if you already have a CDS files)
```
Step 1. Need to have a CDS file with predicted genes 
Step 2. Index the reference (CDS)
Step 3. Run salmon quant
```
#Step 2
 salmon index --transcripts cds_from_genomic.fna --index LCC4_index --threads 20


Step 3: Actually running salmon (for single end reads)
# Set the path to the directory containing the FASTQ files
fastq_dir="/scratch/shrinivas/Symbiodinium_Meta_Analysis/Cleaned_Reads/Gierz_2017/cleaned_reads/"

# Set the index path
index_path="/scratch/shrinivas/Symbiodinium_Meta_Analysis/Reference_Transcriptomes/Fugacium_kawagutti/F_kawagutii_index"

# Set the output directory
output_dir="/scratch/shrinivas/Symbiodinium_Meta_Analysis/Cleaned_Reads/Gierz_2017/salmon_output/"

# Loop over each FASTQ file in the directory
for fastq_file in "$fastq_dir"/*.fastq; do
    # Extract the filename without extension
    filename=$(basename "$fastq_file" .fastq)
    
    # Run Salmon quantification for each FASTQ file
    salmon quant --index "$index_path" \
                 --libType A \
                 --unmatedReads "$fastq_file" \
                 -p 30 \
                 --output "$output_dir/${filename}_salmon_quant"
done

# for paired end reads 
# Set the path to the directory containing the FASTQ files
fastq_dir="/scratch/shrinivas/Symbiodinium_Meta_Analysis/Cleaned_Reads/Lima_2020"

# Set the index path
index_path="/scratch/shrinivas/Symbiodinium_Meta_Analysis/Reference_Transcriptomes/Salmon_Indexes/Symbiodinium_microadriaticum_Index"

# Set the output directory
output_dir="/scratch/shrinivas/Symbiodinium_Meta_Analysis/Cleaned_Reads/Lima_2020/Lima_2020_salmon_output"

# Loop over each paired FASTQ file in the directory
for fastq_file1 in "$fastq_dir"/*_1.fastq; do
    # Check if paired-end file exists
    if [ ! -f "$fastq_file1" ]; then
        continue  # Skip if the first mate is missing
    fi
    
    # Extract the filename without extension
    filename=$(basename "$fastq_file1" _1.fastq)
    
    # Get the second mate file
    fastq_file2="$fastq_dir/${filename}_2.fastq"
    
    # Check if the second mate file exists
    if [ ! -f "$fastq_file2" ]; then
        continue  # Skip if the second mate is missing
    fi
    
    # Run Salmon quantification for each paired-end files
    salmon quant --index "$index_path" \
                 --libType A \
                 --mates1 "$fastq_file1" \
                 --mates2 "$fastq_file2" \
                 -p 30 \
                 --output "$output_dir/${filename}_salmon_quant"
done

# The main counts are gonna be in a file called quant.sf
# Rename the files and move to somewhere better 
#!/bin/bash

# Set the directory containing the directories
parent_directory="/path/to/parent_directory"

# Loop over each directory
for directory in "$parent_directory"/*/; do
    # Extract the directory name
    dir_name=$(basename "$directory")
    
    # Check if quant.sf file exists in the directory
    if [ -f "$directory/quant.sf" ]; then
        # Rename quant.sf to directory_name.tsv
        mv "$directory/quant.sf" "$directory/$dir_name.tsv"
        echo "Renamed quant.sf to $dir_name.tsv in $directory"
    else
        echo "quant.sf not found in $directory"
    fi
done

# Script to move 



#!/bin/bash

# Define the target directory to move .sf files
target_dir="new_directory"

# Create the target directory if it doesn't exist
mkdir -p "$target_dir"

# Iterate over directories starting with "SRR"
for dir in SRR*/; do
    # Check if directory contains .sf files
    if [ -n "$(find "$dir" -maxdepth 1 -type f -name '*.sf')" ]; then
        # Move all .sf files to the target directory
        mv "$dir"*.sf "$target_dir/"
        echo "Moved .sf files from $dir to $target_dir/"
    else
        echo "No .sf files found in $dir"
    fi
done

# change the column names so its easier to filter 
find /scratch/shrinivas/Sargassum/Salmon_output/LCC6/quant_files -type f -name '*.tsv' -exec sh -c 'filename=$(basename "{}" .tsv); awk -v filename="$filename" "BEGIN {FS=OFS=\"\t\"} NR==1 {gsub(/NumReads/, \"NumReads_\"filename); gsub(/TPM/, \"TPM_\"filename)} {print}" "{}" > "$filename"_fixed_headers.tsv' \;

# Extract the counts table my doing simple manipulations
# Note this part of the script is in python 
```
conda activate py2 
import pandas as pd
import os

# Directory containing TSV files
directory = "/scratch/shrinivas/Symbiodinium_Meta_Analysis/Counts_tables/Gierz_2017/fixed_headers/2_col"

# Get list of TSV files
tsv_files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith('.tsv')]

# Read all TSV files into a list of dataframes
dfs = [pd.read_csv(file, sep='\t') for file in tsv_files]

# Merge dataframes based on the "Name" column
merged_df = pd.merge(dfs[0], dfs[1], on='Name')
for df in dfs[2:]:
    merged_df = pd.merge(merged_df, df, on='Name')

# Assuming 'df' is the dataframe you already have and 'stats' is the dataframe you want to merge with
# merged_df = pd.merge(df, stats, on='Calc.MW')
# You can replace 'df' and 'stats' with the appropriate dataframes in your code

# Display the merged dataframe
print(merged_df)

output_file = "merged_output.csv"
merged_df.to_csv(output_file, index=False)

print(f"Merged data saved to '{output_file}'")
```
# Also run Eggnogg-mapper against the .pep files 

/home/timothy/programs/eggnog-mapper-2.1.6/emapper.py -i <input> --cpu 48 --itype proteins -o <output> --data_dir /scratch/timothy/databases/eggnog-mapper-rel20211209

# PHASE 2: ASSESSING DATA
```
use the obtained data and the metadata to rename some of the columns so that they make sense. 
```
