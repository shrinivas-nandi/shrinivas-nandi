#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
#source ~/scripts/script_setup.sh
s#et +eu; conda activate py27; set -eu

export PATH="/home/timothy/miniconda3/envs/py27/bin:$PATH"
export PERL5LIB=/home/timothy/miniconda3/envs/py27/lib/5.26.2/x86_64-linux-thread-multi:/home/timothy/miniconda3/envs/py27/lib/5.26.2
export PATH="/home/timothy/programs/RepeatAnalysis/RepeatModeler-2.0.1:$PATH"

PA=10

#### Start Script
while read REF;
do
  run_cmd "BuildDatabase -name $REF.RMDB $REF 1> ${REF}.BuildDatabase.log 2>&1"
  run_cmd "RepeatModeler -database $REF.RMDB -pa $PA -LTRStruct 1> ${REF}.RepeatModeler.log 2>&1"
done <bin.21.fa

/home/timothy/programs/RepeatAnalysis/RepeatModeler-2.0.1/RepeatModeler -database bin.21.fa.RMDB -pa 5 -LTRStruct 1> bin.21.fa.RepeatModeler.log


/home/timothy/programs/RepeatAnalysis/RepeatModeler-2.0.1/BuildDatabase -name bin.21.fa.RMDB bin.21.fa 1> bin.21.fa.BuildDatabase.log

cp -r /scratch/timothy/projects/0045_Paulinella_MG_2022/03_Analysis/2023-08-16/01_Trimming//PaulHet_RNA* /scratch/shrinivas/Paulinella_Metagenome/bin21/Paul_Het_Transcript_Reads
