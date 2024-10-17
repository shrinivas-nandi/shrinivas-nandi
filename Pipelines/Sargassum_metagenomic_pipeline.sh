Sargassum_metagenomic_pipeline.sh

# QC was done using the naive_atlas pipeline 
naive_atlas run all 

# Assembly (done using megahit on amarel)
megahit -1 /scratch/sn809/Sargassum/reads/QC.errorcorr_R1.fastq.gz -2 /scratch/sn809/Sargassum/reads/QC.errorcorr_R2.fastq.gz -t 40 -o megahit_result --presets meta-large --memory 950

# Binning done using metabat2
metabat2 -i /scratch/sn809/Sargassum/reads/metabat2_wd/final.contigs.fa -o /scratch/sn809/Sargassum/reads/metabat2_wd/metabat2_contig_output.fasta --minContig 1500 --unbinned unbinned.fasta -t 40  


# This pipeline just contains basic commands on how certain procedures were done. And the code that was used to analyze it.  

# from generated files compute basic numbers 
grep -c ">" *.fa > bin_contig_data.csv

# remove bins with too few contigs (<= 10) [optional]

# predict proteins for each (PRODIGAL)

#!/bin/bash

# Define directories
input_dir="/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/downstream/MAGs_of_interest/Up_DAH/genomes"  # Modify this to your input directory path
output_dir="/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/downstream/MAGs_of_interest/Up_DAH/genomes/prodigal_output"  # Modify this to your desired output directory

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all .fa files in the directory
for fa_file in "$input_dir"/*.fasta; do
    # Extract the base filename without extension
    base_name=$(basename "$fa_file" .fasta)

    # Define output file names based on the base name
    output_gbk="$output_dir/${base_name}_genes.gbk"
    output_protein="$output_dir/${base_name}_protein.faa"

    # Run Prodigal
    prodigal -i "$fa_file" -o "$output_gbk" -a "$output_protein"
done

echo "Prodigal processing complete."



# now compute basic stats for each protein 
grep -c ">" *.fa > bin_contig_data.csv

#### Use this code to import files from server to your own computer (an exmaple)
scp -r shrinivas@coral.rutgers.edu:/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/bins_prot.csv  /Users/shrinivas/Desktop




#### extract the longest sequence from each file and move into one central fasta file ####################


#!/bin/bash

# Define directories
input_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files"  # Modify this to your input directory path
output_file="long_seq.fa"  # Output file for longest sequences

# Clear the output file if it exists
> "$output_file"

# Loop through all .fa files in the directory
for fa_file in "$input_dir"/*.faa; do
    # Extract the base filename without extension
    base_name=$(basename "$fa_file" .faa)
    
    # Variables to store the longest sequence and header
    longest_seq=""
    longest_header=""
    longest_length=0

    current_seq=""
    current_header=""

    # Read the .fa file line by line
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            # When a new header is found, check the length of the previous sequence
            if [[ -n $current_seq ]]; then
                seq_length=${#current_seq}
                if [[ $seq_length -gt $longest_length ]]; then
                    longest_seq=$current_seq
                    longest_header=$current_header
                    longest_length=$seq_length
                fi
            fi
            # Start new sequence
            current_header=$line
            current_seq=""
        else
            # Append the sequence lines (removing spaces)
            current_seq+=$line
        fi
    done < "$fa_file"

    # After the last sequence in the file
    if [[ -n $current_seq ]]; then
        seq_length=${#current_seq}
        if [[ $seq_length -gt $longest_length ]]; then
            longest_seq=$current_seq
            longest_header=$current_header
            longest_length=$seq_length
        fi
    fi

    # Write the longest sequence to the output file with the modified header
    echo ">$base_name|$longest_header" >> "$output_file"
    echo "$longest_seq" >> "$output_file"

done

echo "Longest sequences have been written to $output_file."


###### Blastp analysis for the top hit (https://github.com/bbuchfink/diamond)

/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond blastp --ultra-sensitive --iterate --max-target-seqs 0 --evalue 0.00001 --query long_seq.fa  --db /scratch/databases/ncbi/2022_07/nr.dmnd --out blast_output_long_seq --threads 50 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms 
# let it run as a screen


### Run CheckM for identifyin MAG lineage (https://github.com/Ecogenomics/CheckM)

checkm lineage_wf /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/original_files/*.faa /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_file/checkm_output --genes --threads 50 -x .faa


## Run coverm for relative abundance of each MAG 

# samtools 
#!/bin/bash

# Directory where your .bam files are stored
INPUT_DIR="/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/naive_atlas/genomes/alignments/bams"
# Samtools binary path
SAMTOOLS="/home/timothy/programs/samtools-1.11/bin/samtools"
# Number of threads
THREADS=50

# Loop through all .bam files in the directory
for BAM_FILE in ${INPUT_DIR}/*.bam
do
    # Get the base name without .bam extension
    BASE_NAME=$(basename ${BAM_FILE} .bam)
    
    # Define the output sorted file name
    OUTPUT_FILE="sorted_${BASE_NAME}.bam"
    
    # Run samtools sort
    ${SAMTOOLS} sort ${BAM_FILE} -o ${OUTPUT_FILE} --threads ${THREADS}
    
    echo "Sorted ${BAM_FILE} to ${OUTPUT_FILE}"
done


/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/naive_atlas/genomes/alignments/bams


coverm genome --bam-files /scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/naive_atlas/genomes/alignments/bams/sorted_N39.bam --genome-fasta-directory /scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/naive_atlas/genomes/genomes --threads 50 -x .fasta -o /scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/naive_atlas/coverm/N39_coverm_output.csv


#### Run kraken for originial reads 

kraken2  --db /scratch/shrinivas/Databases/kraken_db_standard/ --threads 60  --quick --use-names --paired <>_R1.fastq.gz <>_R2.fastq.gz --report DS_kraken_report.txt --output DS_kraken_out


#!/bin/bash

# Directory where your fastq.gz files are stored
INPUT_DIR="/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/naive_atlas/QC/reads"
# Database location
DB="/scratch/shrinivas/Databases/kraken_db_standard/"
# Number of threads
THREADS=60

# Loop through all *_R1.fastq.gz files in the directory
for R1_FILE in ${INPUT_DIR}/*_R1.fastq.gz
do
    # Get the base name without _R1.fastq.gz
    BASE_NAME=$(basename ${R1_FILE} _R1.fastq.gz)
    
    # Define corresponding R2 file
    R2_FILE=${INPUT_DIR}/${BASE_NAME}_R2.fastq.gz
    
    # Output files
    REPORT="${BASE_NAME}_kraken_report.txt"
    OUTPUT="${BASE_NAME}_kraken_output.txt"
    
    # Run Kraken2 with the paired files
    kraken2 --db ${DB} --threads ${THREADS} --quick --use-names \
            --paired ${R1_FILE} ${R2_FILE} \
            --report ${REPORT} --output ${OUTPUT}
    
    echo "Processed ${BASE_NAME}"
done


#!/bin/bash

# Define paths
input_file="/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/downstream/MAGs_of_interest/Up_DAH/MAGs_of_interest.txt"
source_dir="/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/naive_atlas/genomes/annotations/genes/"  # Ensure no leading space
destination_dir="/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/downstream/MAGs_of_interest/Up_DAH/metaeuk_gene_predictions/"

# Ensure the destination directory exists
mkdir -p "$destination_dir"

# Loop through each line of the input file (filenames without extensions)
while IFS= read -r filename
do
    # Check and copy the .fna file if it exists
    if [ -f "${source_dir}${filename}.fna" ]; then
        echo "Copying ${filename}.fna to $destination_dir"
        cp "${source_dir}${filename}.fna" "$destination_dir"
    else
        echo "${filename}.fna not found."
    fi
    
    # Check and copy the .FAA file if it exists
    if [ -f "${source_dir}${filename}.faa" ]; then
        echo "Copying ${filename}.faa to $destination_dir"
        cp "${source_dir}${filename}.faa" "$destination_dir"
    else
        echo "${filename}.faa not found."
    fi

done < "$input_file"

echo "File copy process completed."


# for eukaryotes pair with metaeuk for protein annotation

metaeuk easy-predict cd/referenceDB predsResults tempFolder


### Use mmseqs2 to annotate these proteins using nr database 
# example code 
mmseqs easy-taxonomy /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/original_files/*.faa /scratch/shrinivas/Databases/mmseq_nr_db/mmseq_nr /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/original_files/mmseq_output /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/original_files/tmp --threads 50


#### Functional annotation of bacteria####
# option 1: eggnog mapper 

#!/bin/bash

# Define the command template and input/output directories
input_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/prokaryotic_mags"
output_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/prokaryotic_mags/eggnon-mapper_annotation"
eggnog_mapper="/home/timothy/programs/eggnog-mapper-2.1.6/emapper.py"
data_dir="/scratch/timothy/databases/eggnog-mapper-rel20211209"

# Ensure output directory exists
mkdir -p "$output_dir"

# Loop through all .faa files in the input directory
for file in "$input_dir"/*.faa; do
    # Extract the base name of the file (without directory and extension)
    base_name=$(basename "$file" .faa)
    
    # Define the output file name
    output_file="$output_dir/${base_name}_eggnog"

    # Run the EggNOG-mapper command
    $eggnog_mapper -i "$file" --cpu 48 --itype proteins -o "$output_file" --data_dir "$data_dir"
    
    echo "Processed $file -> $output_file"
done



# option 2: blast everything
#!/bin/bash

# Define the path to Diamond, database, input, output directories, and log file
diamond_path="/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond"
db_path="/scratch/databases/ncbi/2022_07/nr.dmnd"
output_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/prokaryotic_mags/blast_output"
input_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/prokaryotic_mags"
log_file="$output_dir/diamond_run.log"

# Ensure output directory exists
mkdir -p "$output_dir"

# Initialize log file
echo "Diamond BLAST run started on $(date)" > "$log_file"

# Loop through all .faa files in the input directory
for file in "$input_dir"/*.faa; do
    # Extract the base name of the file (without directory and extension)
    base_name=$(basename "$file" .faa)
    
    # Define the unique output file name
    output_file="$output_dir/${base_name}_diamond_output.txt"
    
    # Log the file being processed
    echo "Processing $file -> $output_file" >> "$log_file"
    
    # Run the Diamond command
    "$diamond_path" blastp --ultra-sensitive --iterate --max-target-seqs 1 --evalue 0.00001 \
    --query "$file" --db "$db_path" --out "$output_file" --threads 80 --outfmt 6 \
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms
    
    # Log completion of the file processing
    echo "Finished processing $file at $(date)" >> "$log_file"
done

# Log completion of the entire run
echo "Diamond BLAST run completed on $(date)" >> "$log_file"



# things to do 
1. Internal blast against sargassum genome-fasta-directory
2. Coverm analysis of bin coverage 
3. Annotation and lyase search
4. 

# Internal Blast 
 makeblastdb
/home/timothy/programs/ncbi-blast-2.13.0+/bin/blastn -db /scratch/shrinivas/Sargassum/reference_genomes/GCSG_db/GCSG_db -out final_contigs_blast_output.tsv -evalue 0.00001 -query final.contigs.fa -max_hsps 1 -num_threads 70 -outfmt 6 qseqid staxids bitscore

#!/bin/bash

# Define the directory containing the input files
input_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/metabat2_wd/euk_mags"

# Define the output directory
output_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/blast_output_nucleotide_euks"

# Define the BLAST database
blast_db="/scratch/shrinivas/Sargassum/reference_genomes/GCSG_db/GCSG_db"

# Define BLAST command options
evalue="0.00001"
threads="50"
outfmt="6"

# Loop through each file in the input directory
for input_file in "$input_dir"/*; do
  # Get the base name of the file (without the path and extension)
  base_name=$(basename "$input_file")

  # Construct the output file path
  output_file="$output_dir/${base_name}_output.tsv"

  # Run BLASTn command
  /home/timothy/programs/ncbi-blast-2.13.0+/bin/blastn \
    -db "$blast_db" \
    -out "$output_file" \
    -evalue "$evalue" \
    -query "$input_file" \
    -max_hsps 1 \
    -num_threads "$threads" \
    -outfmt "$outfmt"

  echo "Processed: $input_file -> $output_file"
done


awk -F'\t' '$10 != "NA" {print $1}' quality_summary.tsv > MAG1646_headers.txt

 /home/timothy/programs/seqkit_v2.3.1/seqkit grep -f /scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/viral_check/checkv_out_MAG1646/MAG1646_headers.txt MAG1646.fasta > MAG01646_put_viral_seq.fasta

prodigal -i MAG01646_put_viral_seq.fasta -o MAG01646_put_viral_genes.gbk -a prot_MAG01646_put_viral_seq.faa

/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond blastp --ultra-sensitive --iterate --max-target-seqs 0 --evalue 0.00001 --query prot_MAG01646_put_viral_seq.faa --db /scratch/databases/ncbi/2022_07/nr.dmnd --out prot_MAG1646_blast_output.tsv --threads 70 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms # the actual blast





# Bam files to make stuff (bowtie and samtools)

bowtie2-build contigs.fa contig_index --threads 20

bowtie2 -x contig_index -1 sample_R1.fastq -2 sample_R2.fastq -S output.sam -p 50

bowtie2 -x /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/bowtie/index_bowtie/index_contigs.fa -1 /scratch/shrinivas/Sargassum/Dominican_Republic/cleaned_reads/S1.errorcorr_R1.fastq.gz -2 /scratch/shrinivas/Sargassum/Dominican_Republic/cleaned_reads/S1.errorcorr_R2.fastq.gz -S final_contigs.sam -p 50


samtools view -S -b aligned_output.sam > aligned_output.bam

samtools sort aligned_output.bam -o aligned_output_sorted.bam

# blobtools
/home/shrinivas/Programs/blobtools/blobtools create -i /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/final.contigs.fa -b /scratch/shrinivas/Sargassum/Dominican_Republic/alignment_files/sorted_final_contigs_aligned.bam -t /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/blobtools/blast.out -o /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/blobtools/test && \
/home/shrinivas/Programs/blobtools/blobtools view -i /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/blobtools/test/test.blobDB.json && \
/home/shrinivas/Programs/blobtools/blobtools plot -i /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/blobtools/test/test.blobDB.json






#checkv complete data
/scratch/shrinivas/Coral_Projects/Dominican_Republic/YBD/viral_check/checkv_genome/checkv_out_MAG1646/blast_output
MAG0768

/home/shrinivas/Programs/blobtools/blobtools create -i /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/final.contigs.fa -b /scratch/shrinivas/Sargassum/Dominican_Republic/alignment_files/sorted_bam.bam --db /home/shrinivas/Programs/blobtools/data/nodesDB.txt -o blob_db_sargassum_only -t /scratch/shrinivas/Sargassum/Dominican_Republic/alignment_files/blobtools_files/sargassum_blobtools_blast.tsv


/home/timothy/programs/ncbi-blast-2.13.0+/bin/blastn -db /scratch/databases/ncbi/2022_07/nt -out contigs_blast_blobtools_output.tsv -evalue 0.00001 -query final.contigs.fa -max_hsps 1 -num_threads 50 -outfmt 6 qseqid staxids bitscore #taxid for blobtools


#########Tasks tbd before 5/10##########
'''
1. Blobtools using sargassum db
2. Bolbtools using regular blast against
3. Pull sargassum reads using tiara and blobtools and ID sargassum bins
4. General lyase and sulfatase across bins (eggnogg)
5. CaZy blast 
6. Heatmap plot
7. GTDB search
'''


#####3. CaZy blast ##
#!/bin/bash

# Paths to Diamond, database, input, output directories, and log file
diamond_path="/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond"
db_path="/scratch/shrinivas/Databases/CaZy_db/CaZy_db.dmnd"
output_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/prokaryotic_mags/cazy_blast"
input_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/prokaryotic_mags"
log_file="$output_dir/diamond_run.log"

# Ensure output directory exists
mkdir -p "$output_dir"

# Initialize log file
echo "Diamond BLAST run started on $(date)" > "$log_file"

# Loop through all .faa files and run them in parallel
find "$input_dir" -name "*.faa" | parallel -j 8 '
    file={}
    base_name=$(basename "$file" .faa)
    output_file="'$output_dir'/${base_name}_diamond_output.txt"
    
    echo "Processing $file -> $output_file" >> "'$log_file'"

    "'$diamond_path'" blastp --ultra-sensitive --iterate --max-target-seqs 1 --evalue 0.00001 \
        --query "$file" --db "'$db_path'" --out "$output_file" --threads 80 --outfmt 6 \
        qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms

    echo "Finished processing $file at $(date)" >> "'$log_file'"
'

# Log completion of the entire run
echo "Diamond BLAST run completed on $(date)" >> "$log_file"

echo "Diamond BLAST run completed on $(date)" >> "$log_file"


# SulfAtlas Blast

#!/bin/bash

# Paths to Diamond, database, input, output directories, and log file


diamond_path="/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond"
db_path="/scratch/shrinivas/Databases/sulfatlas_db/sulfatlas_db.dmnd"
output_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/prokaryotic_mags/sulfatlas_blast"
input_dir="/scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/prodigal_output/protein_files/prokaryotic_mags"
log_file="$output_dir/diamond_run.log"

# Ensure output directory exists
mkdir -p "$output_dir"
cd
# Initialize log file
echo "Diamond BLAST run started on $(date)" > "$log_file"

# Loop through all .faa files and run them in parallel
find "$input_dir" -name "*.faa" | parallel -j 8 '
    file={}
    base_name=$(basename "$file" .faa)
    output_file="'$output_dir'/${base_name}_diamond_output.txt"
    
    echo "Processing $file -> $output_file" >> "'$log_file'"

    "'$diamond_path'" blastp --ultra-sensitive --iterate --max-target-seqs 1 --evalue 0.00001 \
        --query "$file" --db "'$db_path'" --out "$output_file" --threads 80 --outfmt 6 \
        qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms

    echo "Finished processing $file at $(date)" >> "'$log_file'"
'

# Log completion of the entire run
echo "Diamond BLAST run completed on $(date)" >> "$log_file"




/home/shrinivas/Programs/blobtools/blobtools create -i /scratch/shrinivas/Sargassum/Dominican_Republic/S1_bins/final.contigs.fa -b /scratch/shrinivas/Sargassum/Dominican_Republic/alignment_files/sorted_bam.bam \
--db /home/shrinivas/Programs/blobtools/data/nodesDB.txt -o blob_db_sargassum -t /scratch/shrinivas/Sargassum/Dominican_Republic/alignment_files/blobtools_files/contigs_blast_blobtools_output.tsv




#internal blast against sargassum for blobtools

/home/timothy/programs/ncbi-blast-2.13.0+/bin/blastn -db /scratch/shrinivas/Sargassum/reference_genomes/GCSG_db/GCSG_db -out final_contigs_blast_output.tsv -evalue 0.00001 -query final.contigs.fa -max_hsps 1 -num_threads 70 -outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen



##### formatting blast outputs for counts #######

for file in *.txt; do
    awk 'BEGIN {OFS="\t"; print "", "Column2", "CaZyme"} { gsub(/\|/, "\t"); print }' "$file" > "${file%.txt}.tsv"
done


## generating a summary table (example only for sulfatase. Cazyme one was generated in a similar manner)
# Run as a pandas script
import os
import pandas as pd

# Define the directory containing your TSV files
input_dir = '/Users/shrinivas/Desktop/Sargassum/sulfatase_tsv/'
output_dir = '/Users/shrinivas/Desktop/Sargassum/sulfatase_tsv/Sult_counts/' # Define where to save the output files

# Loop through each file, process it, and save the result as a new file
for tsv_file in tsv_files:
    file_path = os.path.join(input_dir, tsv_file)
    
    try:
        # Read the TSV file, assuming it has no headers and process lines with errors
        df = pd.read_csv(file_path, sep='\t', header=None)
        
        # Split the second column (index 1) based on the "|" character
        split_columns = df.iloc[:, 1].str.split('|', expand=True)
        
        # Extract the 'sult' from the second part of the split (index 1)
        df['sult'] = split_columns[1].str.extract(r'_([^_]+_[^_]+)$')
        
        # Generate counts of unique values in the third column (assuming index 2 contains CaZyme data)
        cazyme_column = df.iloc[:, 14].astype(str)  # Ensure the third column is a string
        unique_counts = cazyme_column.value_counts()

        # Convert to DataFrame for better visualization
        unique_counts_df = unique_counts.reset_index()
        # Rename the columns: 'CaZyme' and use the file name for the count column
        unique_counts_df.columns = ['CaZyme', tsv_file]

        # Save the counts file with the appropriate name
        output_file_path = os.path.join(output_dir, f"{os.path.splitext(tsv_file)[0]}_counts.csv")
        unique_counts_df.to_csv(output_file_path, index=False)
        print(f"Processed and saved counts file: {output_file_path}")
    
    except pd.errors.ParserError as e:
        print(f"Error processing file {tsv_file}: {e}")
    except Exception as e:
        print(f"An error occurred with file {tsv_file}: {e}")


# merge togehter 
# merge all the csv on CaZyme
# Define the directory containing your CSV files
csv_dir = '/Users/shrinivas/Desktop/Sargassum/sulfatase_tsv/Sult_counts/'
output_file = '/Users/shrinivas/Desktop/Sargassum/sulfatase_counts.csv'

# List all CSV files in the directory
csv_files = [f for f in os.listdir(csv_dir) if f.endswith('.csv')]

# Initialize an empty DataFrame to store the merged results
merged_df = None

# Loop through each CSV file and merge them on the 'CaZyme' column
for csv_file in csv_files:
    file_path = os.path.join(csv_dir, csv_file)
    
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Merge on 'CaZyme'. If merged_df is None, this is the first file, so we start with it.
    if merged_df is None:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on='CaZyme', how='outer')

# Save the final merged DataFrame to a CSV file
merged_df.to_csv(output_file, index=False)

print(f"All CSV files merged and saved to: {output_file}")





#####

awk 'BEGIN {FS=OFS="\t"} {print $0, 3015}' final_contigs_blast_output.tsv > output_file.tsv

cut -f1,15,12 output_file.tsv > sargassum_blobtools_blast.tsv



###### Jellyfish to estimate genome size from reads #######

/home/timothy/programs/jellyfish-2.3.0/bin/jellyfish count --mer-len 21 -s 80 --threads 70 --output 21_mer_out -C /scratch/shrinivas/Sargassum/Dominican_Republic/cleaned_reads/S1.errorcorr_R1.fastq /scratch/shrinivas/Sargassum/Dominican_Republic/cleaned_reads/S1.errorcorr_R2.fastq





