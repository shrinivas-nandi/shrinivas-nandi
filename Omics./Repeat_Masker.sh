#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1


export PATH="$PATH:/home/timothy/programs/RepeatAnalysis/RepeatMasker-4.1.2-p1"
export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
export PATH="$PATH:/home/timothy/programs/bbmap"

UTILS="/home/timothy/programs/RepeatAnalysis/RepeatMasker-4.1.2-p1/util"
PA=12

#### Start Script
while read REF;
do
  echo -e "\n\n\n\n\n"
  log "## bin.21.fa"
  LIB="${REF}.custom_repeat_library.fa"
  run_cmd "cat ${REF}.RMDB-families.fa rmlib.fa > ${LIB}"
  run_cmd "RepeatMasker -no_is -a -x -gff -pa ${PA} -lib ${LIB} ${REF} 1>${REF}.RepeatMasker.log 2>&1"
  run_cmd "bedtools maskfasta -fi ${REF} -fo ${REF}.softmasked -bed ${REF}.out.gff -soft"
done


/home/timothy/programs/bedtools-2.29.2/bin/bedtools maskfasta -fi bin.21.fa.masked -fo bin.21.fa.softmasked -bed bin.21.fa.out.gff  -soft
