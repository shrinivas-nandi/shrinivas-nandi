STAR 
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