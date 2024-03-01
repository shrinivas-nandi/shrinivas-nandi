hisat2.sh 
# build a reference
hisat2-build -f Symbiodinium_kawagutii.assembly.935Mb.fa ./Alt_ref

1. Align against reference genome 
$F=/scratch/shrinivas/Microbiome_Chille/cat/
hisat2 -p 8 --rf --dta -q -x Reference_Genomes/Microbiome_Ref -1 cleaned_reads/Microbio9_S110_R1_001.clean.fastq.gz -2 cleaned_reads/Microbio9_S110_R2_001.clean.fastq.gz -S hisat2/Microbio9_S110.sam
