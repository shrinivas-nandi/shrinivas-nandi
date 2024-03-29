#!/usr/bin/env bash

# Metagenome Analysis

This code was developed and utilzied for the work I did on the P.ovalis data in 2023 in the Bhattacharya Lab.
I will eventually add Timothy Stephens to edit this document to add snippets of code.

# Phase I: Generating the Bins

## Step 1: QC of raw reads obtained from sequencing 
Tools used 
- Fastqc 
- Multiqc
Create a directory with the raw reads only 

``` 
fastqc . 
# subsequently ran multiqc 
multiqc .
```

## Step 2: Cleaning the reads 
Tools used 
- Fastp 

```
sh -c 'for file in "Microbio1_S103" "Microbio3_S105" "Microbio4_S106" "Microbio5_S107" "Microbio6_S108" "Microbio8_S109" "Microbio9_S110"

do
fastp --in1 ${file}_R1_001.fastq.gz --in2 ${file}_R2_001.fastq.gz\
 --out1 /scratch/shrinivas/Microbiome_Chille/cat/cleaned_reads/${file}_R1_001.clean.fastq.gz --out2 /scratch/shrinivas/Microbiome_Chille/cat/cleaned_reads/${file}_R2_001.clean.fastq.gz\
 --failed_out /scratch/shrinivas/Microbiome_Chille/cat/cleaned_reads/${file}_failed.txt \
--qualified_quality_phred 20 --unqualified_percent_limit 10 --length_required 100\
 detect_adapter_for_pe --cut_right cut_right_window_size 5 cut_right_mean_quality 20
done'
```

## Step 3: Run qc on the cleaned reads
``` 
fastqc . 
# subsequently ran multiqc 
multiqc .
```

## Step 4: Align reads
Tools used 
- Bowtie: Better for Short reads. Commonly used for microbiome studies
or 
- STAR: Good for transcriptomics and obtaining a counts table. But I also used this tool for aligning stuff in the MAG paper so I'm not sure
Note: STAR output is also required for running BRAKER in the future.

**Question**: What is the difference between Hisat2 and STAR? Which aligner is better? 

```
# Bowtie2
/home/timothy/programs/bowtie2-2.4.4-linux-x86_64/bowtie2 -x Paul_Het_index -1 PaulHet_DNA_Normal_L1_trimmed_R1.fastq.gz -2 PaulHet_DNA_Normal_L1_trimmed_R2.fastq.gz -p 50  > Paul_Het_Aligned.sam 2> Paul_Het_Aligned.stats

#STAR
/home/timothy/programs/STAR-2.7.8a/bin/Linux_x86_64_static/STAR --genomeDir /scratch/shrinivas/Paulinella_Metagenome/bin21/Paul_Het_Transcript_Reads --genomeFastaFiles bin.21.fa --readFilesIn "PaulHet_RNA_L1_trimmed_R1.fastq.gz"  "PaulHet_RNA_L1_trimmed_R2.fastq.gz" --outFileNamePrefix output_directory
```

## Step 5: Splitting the reads into Bins 
**Notes**: Tim what software did you use for binning 


## Step 6: Obtain basic stats on the bins 
- bbmaps

```
# need to input code for this
```


## Step 7: Coverm to see what % of assembly to mapped to a bin 
```
need to add the input code
```

# Phase II: Working on the bins 

## Step 1: Find out what kingdom is in each bin 
Tools used 
Busco 
```
# checkm 
checkm lineage_wf /scratch/shrinivas/Paulinella_Metagenome/PaulHet_DNA.spades.scaffolds.fasta.MetaBAT2_bins/bins_fa /scratch/shrinivas/Paulinella_Metagenome/checkM_Attempt2 --reduced_tree -t 30 -x .fa
# busco
/scratch/shrinivas/Paulinella_Metagenome/Busco/busco/bin/busco -i EP00466_Bigelowiella_natans.fasta -o busco_Bigelowiella_natans -m prot -l eukaryota_odb10
```

## Step 2: Predict Open reading frames: Basically split into amino acid codons 
- Tools used: PRODIGAL
```
prodigal  -i bin.10.fa -o output_genes.gbk -a output_protein.faa
```

## Step 3: Run a BLASTp to analyze what is in each bin 
**Take code from Blast file**


## Step 4: Predict Orthogroups with an Outgroup 
```
/home/timothy/programs/OrthoFinder_v2.5.4/orthofinder -f /scratch/shrinivas/Paulinella_Metagenome/Oct25_2023/All_Euk_Tree
```

# PHASE III: Genome Clean Up 
Assuming there is bin of interest from which we want to clean up a genome for publishing and future analysis, like bin21 in P.ovalis experiment. 

## Step 1: Check genome completeness 
- Tools used: BUSCO 
```
/scratch/shrinivas/Paulinella_Metagenome/Busco/busco/bin/busco -i EP00466_Bigelowiella_natans.fasta -o busco_Bigelowiella_natans -m prot -l eukaryota_odb10
```
