Braker.sh 

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

