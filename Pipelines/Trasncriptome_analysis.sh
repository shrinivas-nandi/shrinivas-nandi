hisat2-build -f ../ref/Mcap.genome_assembly.fa ./Mcap_ref


$F=/scratch/shrinivas/Transcript/cleaned_reads/ hisat2 -p 8 --rf --dta -q -x Reference_Genome/Symb_ref -1 SRR1300302.fastq.gz -S hisat2/SRR1300302.sam


hisat2-build -f Symbiodinium_kawagutii.assembly.935Mb.fa.gz ./Alt_ref # always run against a genome 

sh -c 'for file in "SRR17809146" "SRR17809145" "SRR17809144" "SRR17809157" "SRR17809156"
"SRR17809155"
"SRR17809154"
"SRR17809147"
"SRR14297391"
"SRR14297383"
"SRR14297392"
"SRR14297402"
"SRR14297403"
"SRR14297372"
"SRR14297373"
"SRR14297390"
"SRR14297382"
"SRR14297381"
"SRR14297380"
"SRR14297361"
"SRR14297362"
"SRR14297363"
"SRR14297364"
"SRR14297401"

do
fastp --in1 ${file}_1.fastq --in2 ${file}_2.fastq\
 --out1 /scratch/shrinivas/Transcript/Camp_Paper/Breviolum/cleaned_reads/${file}_1.clean.fastq --out2 /scratch/shrinivas/Transcript/Camp_Paper/Breviolum/cleaned_reads/${file}_2.clean.fastq\
 --failed_out /scratch/shrinivas/Transcript/Camp_Paper/Breviolum/cleaned_reads/${file}_failed.txt \
--qualified_quality_phred 20 --unqualified_percent_limit 10 --length_required 75\
 detect_adapter_for_pe --cut_right cut_right_window_size 5 cut_right_mean_quality 20
done'


for file in "SRR17809165" "SRR17809161" "SRR17809160" "SRR17809159" "SRR17809164" \
    "SRR17809158" "SRR17809153" "SRR17809148" "SRR14297370" "SRR14297397" \
    "SRR14297398" "SRR14297399" "SRR14297371" "SRR14297400" "SRR14297388" \
    "SRR14297389" "SRR14297357" "SRR14297385" "SRR14297358" "SRR14297359" \
    "SRR14297356" "SRR14297360" "SRR14297365" "SRR14297384"

do
    fastp --in1 "${file}_1.fastq" --in2 "${file}_2.fastq" \
        --out1 "/scratch/shrinivas/Transcript/Camp_Paper/Durusdinium/cleaned_reads/${file}_1.clean.fastq" \
        --out2 "/scratch/shrinivas/Transcript/Camp_Paper/Durusdinium/cleaned_reads/${file}_2.clean.fastq" \
        --failed_out "/scratch/shrinivas/Transcript/Camp_Paper/Durusdinium/cleaned_reads/${file}_failed.txt" \
        --qualified_quality_phred 20 --unqualified_percent_limit 10 --length_required 75 \
        --detect_adapter_for_pe --cut_right --cut_right_window_size 5 --cut_right_mean_quality 20
done


############### 1_09_2023


sh -c 'for file in "SRR17809146" "SRR17809145" "SRR17809144" "SRR17809157" "SRR17809156"
"SRR17809155"
"SRR17809154"
"SRR17809147"
"SRR14297391"
"SRR14297383"
"SRR14297392"
"SRR14297402"
"SRR14297403"
"SRR14297372"
"SRR14297373"
"SRR14297390"
"SRR14297382"
"SRR14297381"
"SRR14297380"
"SRR14297361"
"SRR14297362"
"SRR14297363"
"SRR14297364"
"SRR14297401"

do
fastp --in1 ${file}_1.fastq --in2 ${file}_2.fastq\
 --out1 /scratch/shrinivas/Transcript/Camp_Paper/Breviolum/cleaned_reads/${file}_1.clean.fastq --out2 /scratch/shrinivas/Transcript/Camp_Paper/Breviolum/cleaned_reads/${file}_2.clean.fastq\
 --failed_out /scratch/shrinivas/Transcript/Camp_Paper/Breviolum/cleaned_reads/${file}_failed.txt \
--qualified_quality_phred 20 --unqualified_percent_limit 10 --length_required 75\
 detect_adapter_for_pe --cut_right cut_right_window_size 5 cut_right_mean_quality 20
done'


### doing salmon 

sh -c 'for file in "SRR17809146" "SRR17809145" "SRR17809144" "SRR17809157" "SRR17809156"
"SRR17809155"
"SRR17809154"
"SRR17809147"
"SRR14297391"
"SRR14297383"
"SRR14297392"
"SRR14297402"
"SRR14297403"
"SRR14297372"
"SRR14297373"
"SRR14297390"
"SRR14297382"
"SRR14297381"
"SRR14297380"
"SRR14297361"
"SRR14297362"
"SRR14297363"
"SRR14297364"
"SRR14297401"

do
salmon quant -i /scratch/shrinivas/Transcript/Camp_Paper/Genomes/salmon_index_Breviolum -l A -1 {file}_1.clean.fastq -2 {file}_2.clean.fastq -o /scratch/shrinivas/Transcript/Camp_Paper/Breviolum/cleaned_reads/salmon_output_directory
done'

# salmon index

salmon index -t Durusdinium_trenchii_CCMP2556_v1.topIsoform.cds.fna -i salmon_index_Durusdinium



# salmon_index for Cladocopium
for file in "SRR17809143" "SRR17809142" "SRR17809163" "SRR17809162" "SRR17809152" \
    "SRR17809151" "SRR17809150" "SRR17809149" "SRR14297393" "SRR14297394" \
    "SRR14297395" "SRR14297396" "SRR14297374" "SRR14297375" "SRR14297376" \
    "SRR14297377" "SRR14297379" "SRR14297378" "SRR14297387" "SRR14297386" \
    "SRR14297366" "SRR14297367" "SRR14297368" "SRR14297369"; do

  # Create a unique output directory for each SRR
  output_dir="/scratch/shrinivas/Transcript/Camp_Paper/Cladocopium/cleaned_reads/salmon_output_directory_revised/$file"
  mkdir -p "$output_dir"
  
  # Run salmon quant command
  salmon quant -i "/scratch/shrinivas/Transcript/Camp_Paper/Genomes/salmon_index" -l A -1 "${file}_1.clean.fastq" -2 "${file}_2.clean.fastq" -o "$output_dir"
done

# for Durusdinium


scp -r ssh shrinivas@coral.rutgers.edu:/scratch/shrinivas/Transcript/Camp_Paper/Durusdinium/cleaned_reads/salmon_output_directory/merged_with_labels.sf /Users/Shrinivas/Desktop

scp -r shrinivas@coral.rutgers.edu:/scratch/shrinivas/Transcript/Camp_Paper/Breviolum/cleaned_reads/salmon_output_directory /Users/Shrinivas/Desktop

/home/timothy/programs/my_interproscan/interproscan-5.53-87.0/interproscan.sh -dp --goterms --cpu 8 --output-file-base SymbC_genes.pep.InterProScan --tempdir SymbCgenes.pep.InterProScan.temp --input Cladocopium_goreaui.PEP.fasta



#!/bin/bash

# Directory containing the subdirectories with .sf files
main_directory=/scratch/shrinivas/Transcript/Camp_Paper/Durusdinium/cleaned_reads/salmon_output_directory

# Output file
output_file="merged.sf"

# Get a list of subdirectory names
subdirectories=$(find "$main_directory" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)

# Initialize the header with directory names
header=$(echo -e "Directory\t$(echo "$subdirectories" | tr '\n' '\t')")

# Add the header to the output file
echo -e "$header" > "$output_file"

# Loop through each subdirectory
for subdir in $subdirectories; do
    # Add the subdirectory name as a column header
    echo -e "$subdir\t$(paste "$main_directory/$subdir"/*.sf | tr '\t' '\n')" >> "$output_file"
done

### Works

hisat2 -p 8 --rf --dta -q -x /scratch/shrinivas/Transcript/Lin_Paper/cleaned_reads/Reference_Genome/Symb_ref.fasta -U /scratch/shrinivas/Transcript/Lin_Paper/cleaned_reads/SRR1300304.fastq.gz -S SRR1300304.sam

## Convert sam files to bam files 

samtools view -bS SRR1300302.sam > SRR1300302.bam


## Conda environment list 
- mutltiqc --> for multiqc, samtools
- Microbiome_Chille --> fastp, fastqc

stringtie sorted_output.bam -G transcripts.gtf -e -B -A transcript_abundances.gtf

awk '$3 == "transcript" { match($0, /transcript_id "([^"]+)"/, tid); match($0, /FPKM "([^"]+)"/, fpkm); match($0, /TPM "([^"]+)"/, tpm); print tid[1], fpkm[1], tpm[1] }' transcript_abundances.gtf > transcript_abundances.txt


### Once its converted into a sam file post indexing 
--> Do hisat2 and align with the reference genome 
hisat2 -p 8 --rf --dta -q -x /scratch/shrinivas/Transcript/Lin_Paper/cleaned_reads/Reference_Genome/Symb_ref.fasta -U /scratch/shrinivas/Transcript/Lin_Paper/cleaned_reads/SRR1300304.fastq.gz -S SRR1300304.sam


### Now convert the sam file to a bam file using samtools
samtools view -bS SRR1300304.sam >> SRR1300304.bam

### sort the bam file 
samtools sort SRR1300304.bam -o sorted_output.bam

### use stringtie to obtain transcripts
stringtie sorted_output.bam -o transcripts.gtf

stringtie -eB -G transcripts.gtf -o transcript_counts.txt sorted_output.bam # maybe 

stringtie sorted_output.bam -G transcripts.gtf -e -B -A transcript_abundances.gtf

# convert that into a text file 
cat transcript_abundances.gtf >> SRR1300304.txt 


cp /scratch/timothy/genome_data/Symbiodinium_Genomes/coral_meta_omics/Durusdinium_trenchii_CCMP2556_v1.topIsoform.pep.faa /scratch/shrinivas/Transcript/Camp_Paper/Genomes

/home/timothy/programs/my_interproscan/interproscan-5.53-87.0/interproscan.sh -dp --goterms --cpu 8 --output-file-base SymbD_genes.pep.InterProScan --tempdir SymbDgenes.pep.InterProScan.temp --input Durusdinium_trenchii_CCMP2556_v1.topIsoform.pep.faa
