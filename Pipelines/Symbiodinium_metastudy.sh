Symbiodinium_metastudy.sh


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

# Clean up the fasta files 
fastp -i SRR4237276.fastq -o SRR4237276_clean.fastq --failed_out SRR4237276_failed_out.txt --qualified_quality_phred 20 --unqualified_percent_limit 10 --thread 16 

# Now once clean. We run SALMON (quasi mapping. doesnt do gene prediction. Good if you already have a CDS files)
```
Step 1. Need to have a CDS file with predicted genes 
Step 2. Index the reference (CDS)
Step 3. Run salmon quant
```
Step 2 (need to finish)
salmon index -h 

Step 3: Actually running salmon 
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

# The main counts are gonna be in a file called quant.sf
# Rename the files and move to somewhere better 

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



