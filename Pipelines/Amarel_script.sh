#!/bin/bash

#SBATCH --partition=main                        # Partition (job queue)

#SBATCH --requeue                               # Re-run job if preempted

#SBATCH --export=NONE                           # Export current environment variables to the launched application

#SBATCH --nodes=1                               # Number of nodes

#SBATCH --ntasks=1                              # Number of tasks (usually = cores) on each node

#SBATCH --cpus-per-task=52                      # Cores per task (>1 if multithread tasks)

#SBATCH --output=%x.slurm_out.%j-%2t.%a         # STDOUT output file (will also contain STDERR if --error is not specified)

#SBATCH --mem=180G                              # Real memory (RAM) per node required (MB)

#SBATCH --time=72:00:00                         # Total run time limit (HH:MM:SS)

#SBATCH --job-name=diamond_blastp_NCBInr        # Replace with your jobname

#SBATCH --array=1-100                           # Specify array range

 
#### Pre-run setup

set -e -o pipefail

source ~/slurm_config_v0.4.sh

prerun_info # Print info

cd $SLURM_SUBMIT_DIR

input_directory="/scratch/sn809/Paul_Het/Protein_Files/split_output"

diamond_bin="/scratch/sn809/Programs/DIAMOND/bin/diamond"

# Set the path to the DIAMOND database
diamond_db="/scratch/sn809/Database/nr.dmnd"



#### Load list of commands/files into array

index=0

while read line ; do

        index=$(($index+1))

        filearray[$index]="$line"

done < files2run.txt

QUERY="${filearray[$SLURM_ARRAY_TASK_ID]}"

echo "Line selected from list of files is ${QUERY}"

OUT="${QUERY}.diamond_blastp_RefSeqComplete.outfmt6"


#### Start Script

run_cmd "diamond blastp --ultra-sensitive --max-target-seqs 500 --evalue 0.001 --query ${QUERY} --db ${diamond_db} --out ${OUT} --compress 1 --threads ${SLURM_CPUS_PER_TASK} --outfmt 6 qseqid sseqid pident length misma

tch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms"

 


 

#### Post-run info

postrun_info # Print info



# for single 
/home/timothy/programs/DIAMOND_v2.0.15/bin/diamond blastp --ultra-sensitive --max-target-seqs 500 --evalue 0.001 --query extracted_sequences.fa --db /scratch/databases/ncbi/2022_07/nr.dmnd --out extracted_sequences_blast.tsv --threads 50 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms
