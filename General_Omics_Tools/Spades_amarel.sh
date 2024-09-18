'''SPADES
spades_use_scaffolds: true # if false use contigs
#Comma-separated list of k-mer sizes to be used (all values must be odd, less than 128 and listed in ascending order).
spades_k: auto
spades_preset: meta # meta, ,normal, rna  single end libraries doesn't work for metaspades
spades_extra: ''
longread_type: none # [none,"pacbio", "nanopore", "sanger", "trusted-contigs", "untrusted-contigs"]
# Preprocessed long reads can be defined in the sample table with 'longreads' , for more info see the spades manual

'''

# downsampling 
/home/timothy/programs/bbmap/bbnorm.sh in=QC.errorcorr_R1.fastq.gz in2=QC.errorcorr_R2.fastq.gz out=S1_QC.fastq.gz outt=excludedreadfile.fastq.gz reads=1.59^10

/home/sn809/Programs/SPAdes-4.0.0-Linux/bin/

in=<>
in2=<>
out=<name_of_file_to_keep>
outt=<excludedreadfile>

reads<1.59^10>

###### actual script########

#!/bin/bash

#SBATCH --partition=mem # Partition (job queue) 

#SBATCH --no-requeue
# Do not re-run job if preempted

#SBATCH --export=NONE 
# Do not export current env to compute nodes

#SBATCH --nodes=1 
# Number of nodes 

#SBATCH --ntasks=1 
# Number of tasks (usually = cores) on each node

#SBATCH --cpus-per-task=40
# Cores per task (>1 if multithread tasks)

#SBATCH --output=%x.slurm_out.%j-%2t
# STDOUT output file (will also contain STDERR if --error is not specified)

#SBATCH --mem=950G 
# Real memory (RAM) per node required (MB) 

#SBATCH --time=72:00:00
# Total run time limit (HH:MM:SS) 

#SBATCH --job-name=bigmem # Replace with your jobname



#### Pre-run setup

cd $SLURM_SUBMIT_DIR



## Run some stuff

# SPAdes command 



echo $SLURM_CPUS_PER_TASK


'''
/home/sn809/Programs/SPAdes-4.0.0-Linux/bin/spades.py -1 /scratch/sn809/Sargassum/down_sampled_reads/QC.errorcorr_downsampled_R1.fastq.gz -2 /scratch/sn809/Sargassum/down_sampled_reads/QC.errorcorr_downsampled_R2.fastq.gz --meta --threads {} -o /scratch/sn809/Sargassum/down_sampled_reads/assembly_output



### not a part of it 
/home/timothy/programs/bbmap/bbnorm.sh in=QC.errorcorr_R1.fastq.gz in2=QC.errorcorr_R2.fastq.gz out=QC.errorcorr_downsampled_R1.fastq.gz out2=QC.errorcorr_downsampled_R2.fastq.gz outt=excludedreadfile.fastq.gz target=30


