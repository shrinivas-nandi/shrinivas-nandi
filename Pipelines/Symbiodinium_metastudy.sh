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



