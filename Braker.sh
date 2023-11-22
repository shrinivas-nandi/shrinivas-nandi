# BRAKER script
#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate BRAKER; set -eu

export GENEMARK_PATH="/home/timothy/programs/BRAKER/gmes_linux_64"
export AUGUSTUS_CONFIG_PATH="/home/timothy/programs/BRAKER/Augustus/config"
export AUGUSTUS_BIN_PATH="/home/timothy/programs/BRAKER/Augustus/bin"
export AUGUSTUS_SCRIPTS_PATH="/home/timothy/programs/BRAKER/Augustus/scripts"
export BAMTOOLS_PATH="/home/timothy/programs/BRAKER/bamtools/bin"
export BLAST_PATH="/home/timothy/programs/ncbi-blast-2.10.1+/bin"
export PROTHINT_PATH="/home/timothy/programs/BRAKER/ProtHint"
export SAMTOOLS_PATH="/home/timothy/programs/samtools-1.11/bin"

export PATH="/home/timothy/programs/BRAKER/BRAKER-2.1.6/scripts:$PATH"
export PATH="/home/timothy/miniconda3/envs/BRAKER/bin:$PATH"
export PERL5LIB="/home/timothy/miniconda3/envs/BRAKER/lib/5.26.2/x86_64-linux-thread-multi:/home/timothy/miniconda3/envs/BRAKER/lib/5.26.2"

NCPUS=48

#### Start Script
SPECIES="Paulinella_micropora_KR01"
REF="bin.21.fa"
REF_SOFT="${REF}.softmasked"
BAM="${REF}.STARAligned.sortedByCoord.out.bam"
run_cmd "braker.pl --genome=${REF_SOFT} --workingdir=${REF}.braker --bam=${BAM} --cores=${NCPUS} --species=${SPECIES} --softmasking --gff3 --min_contig=2000 --skipAllTraining --useexisting"









# STAR to generate a bam file 
# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu
ulimit -n 10000

export PATH="$PATH:/home/timothy/programs/STAR-2.7.8a/bin/Linux_x86_64_static"
export PATH="$PATH:/home/timothy/programs/bbmap"
NCPUS=48

MT_SAMPLES="samples_MetaT.soil.txt"


#### Start Script
while read REF;
do
  SIZE=$(stats.sh ${REF} format=2 | awk '$1=="All" {print $4}')
  SA=$(printf %.0f $(echo "((l(213596) / l(2)) / 2) - 1" | bc -l))
  echo "SA:${SA}"
  run_cmd "STAR --runThreadN ${NCPUS} --runMode genomeGenerate --genomeDir ${REF}.STAR --genomeFastaFiles ${REF} --genomeSAindexNbases ${SA} > ${REF}.star-genomeGenerate.log 2>&1"
  R1=$(awk '{if(NR==1){L="../../../../00_databases/metatranscriptome_PolyA_corr_reads/"$1".PolyA.filtered_1.fastq.gz"}else{L=L",""../../../../00_databases/metatranscriptome_PolyA_corr_reads/"$1".PolyA.filtered_1.fastq.gz"}}END{print L}' "${MT_SAMPLES}")
  R2=$(awk '{if(NR==1){L="../../../../00_databases/metatranscriptome_PolyA_corr_reads/"$1".PolyA.filtered_2.fastq.gz"}else{L=L",""../../../../00_databases/metatranscriptome_PolyA_corr_reads/"$1".PolyA.filtered_2.fastq.gz"}}END{print L}' "${MT_SAMPLES}")
  OUT="${REF}.STAR"
  run_cmd "STAR --limitBAMsortRAM 350000000000 --twopassMode Basic --runThreadN ${NCPUS} --genomeDir ${REF}.STAR --outFileNamePrefix ${OUT} --outSAMtype BAM SortedByCoordinate --readFilesIn ${R1} ${R2} --readFilesCommand zcat"
done < files2run.txt

# note for other cases. Use bbmap stats.sh to see stats of a genome 
# One line code

/home/timothy/programs/STAR-2.7.8a/bin/Linux_x86_64_static/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /scratch/shrinivas/Paulinella_Metagenome/bin21/Paul_Het_Transcript_Reads/bin_index.STAR --genomeFastaFiles bin.21.fa --genomeSAindexNbases 15 > bin.21.star-genomeGenerate.log

/home/timothy/programs/STAR-2.7.8a/bin/Linux_x86_64_static/STAR  --twopassMode Basic --runThreadN 20 --genomeDir bin_index.STAR --outFileNamePrefix Paul_HET_aligned_reads --outSAMtype BAM SortedByCoordinate --readFilesIn PaulHet_RNA_L1_trimmed_R1.fastq.gz,PaulHet_RNA_L2_trimmed_R1.fastq.gz PaulHet_RNA_L1_trimmed_R2.fastq.gz,PaulHet_RNA_L2_trimmed_R2.fastq.gz --readFilesCommand zcat
