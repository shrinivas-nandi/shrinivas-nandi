######################################### Multigene phylogenetic tree prokaryotes ########################################################################

# identify genes
singularity exec \
  -B /scratch/shrinivas/ \
  /scratch/singularity/GTDB-Tk_v2.3.2-rev1.sif \
  gtdbtk identify \
  --genome_dir /scratch/shrinivas/Sargassum/meta-analysis/naive_atlas_output/genomes/prokayotic_genomes \
  -x fa \
  --out_dir /scratch/shrinivas/Sargassum/meta-analysis/naive_atlas_output/genomes/prokayotic_genomes/gtdbtk_identify \
  --cpu 60

# align and skip reference genomes 
singularity exec \
  -B /scratch/shrinivas/ \
  /scratch/singularity/GTDB-Tk_v2.3.2-rev1.sif \
  gtdbtk align \
  --identify_dir /scratch/shrinivas/Sargassum/meta-analysis/naive_atlas_output/genomes/prokayotic_genomes/gtdbtk_identify \
  --out_dir /scratch/shrinivas/Sargassum/meta-analysis/naive_atlas_output/genomes/prokayotic_genomes/gtdbtk_identify/mafft_output \
  --skip_gtdb_refs \
  --cpu 30

# trimaAI
trimal -automated1 -in gtdbtk.bac120.user_msa.fasta.gz -out gtdbtk.bac120.user_msa.aligned.trimmed.fa

#iqtree2
/home/timothy/programs/iqtree-1.6.12-Linux/bin/iqtree -s /scratch/shrinivas/Sargassum/Degradation_experiment/prokaryotic_genomes/phylogeny/align/gtdbtk.bac120.user_msa.aligned.trimmed.fa -m TEST -bb 1000 -nt AUTO -pre /scratch/shrinivas/Sargassum/Degradation_experiment/prokaryotic_genomes/phylogeny/align/iqtree_res

#####################Hmmer CAZyme set####################
"""selected options 
domtblout: domain specific hmms. can have more than one
tblout: sequence specific
"""
hmmscan \
  --cpu 64 \
  --domtblout prokaryotic_MAGs_vs_dbCAN.domtblout \
  --tblout prokaryotic_MAGs_vs_dbCAN.tblout \
  -E 1e-15 --domE 1e-15 \
  /scratch/shrinivas/Databases/CaZy_db/dbCAN-HMMdb-V14.txt \
  /scratch/shrinivas/Sargassum/meta-analysis/predicted_proteins/prokaryotic_proteome/prokaryotic_proteome.faa \
  > allMAGs_vs_dbCAN.log

############################## Assess competeing models ##############################
'''
This would mean decent coverage of the hmmer. dbCan recommends 0.35 so I will stick to that number. 
Basically, if my odmain maps to atelast 35% of the hmmer domain. 
Columns to be used. 
Interpeting some the coordinates
tlen --> lenght of hmm
qlen --> lenght of my full protein 
hmm_cord --> where my hmm aligned to the db hmm 
pos: 1                                            338
     |-----------------------------------------------|
              11                   255
              |---------------------|
                   hmm_from/to

ali_cord --> where on my complete sequence is this domain 

Your sequence (qlen=393):
pos: 1                                                393
     |--------------------------------------------------|
               50                        295
               |-------------------------|
                      ali_from/to
                      
env_cord --> slightly longer region on my sequence 

Your sequence (qlen=393):
pos: 1                                                393
     |--------------------------------------------------|
            45                            315
            |-----------------------------|
                      env_from/to

               50                   295
               |---------------------|
                    ali_from/to (sits inside env)
                    
for our formula I will first calculate coverage 

# lenght of my domain / lenght of model 
(hmm_to - hmm_from + 1) / tlen
'''

awk '!/^#/ {
    domain_coverage = ($17 - $16 + 1) / $3
    print $0, domain_coverage
}' prokaryotic_MAGs_vs_dbCAN.domtblout >> prokaryotic_MAGs_vs_dbCAN.domtblout.tsv
# cut off used 0.60 quite stringent 

################ Compare and look for overlapping domains on the same proteins ######################
import pandas as pd
import sys

# Column positions (0-indexed)
COL_QUERY    = 3   # protein name
COL_IEVALUE  = 12  # i_evalue
COL_ALI_FROM = 17  # ali_from
COL_ALI_TO   = 18  # ali_to
COL_COVERAGE = 22  # domain_coverage


def calc_overlap(ali_from1, ali_to1, ali_from2, ali_to2):
    overlap_start = max(ali_from1, ali_from2)
    overlap_end   = min(ali_to1, ali_to2)
    overlap_len   = max(0, overlap_end - overlap_start + 1)

    shorter = min(ali_to1 - ali_from1 + 1, ali_to2 - ali_from2 + 1)
    return overlap_len / shorter


def resolve_overlaps(df, overlap_threshold=0.30):
    kept = []

    for protein, group in df.groupby(COL_QUERY):
        group = group.reset_index(drop=True)

        if len(group) == 1:
            kept.append(group)
            continue

        to_drop = set()
        rows = list(group.iterrows())

        for i, (idx1, row1) in enumerate(rows):
            if idx1 in to_drop:
                continue
            for j, (idx2, row2) in enumerate(rows):
                if i >= j:
                    continue
                if idx2 in to_drop:
                    continue

                overlap = calc_overlap(
                    row1[COL_ALI_FROM], row1[COL_ALI_TO],
                    row2[COL_ALI_FROM], row2[COL_ALI_TO]
                )

                if overlap >= overlap_threshold:
                    # lower i_evalue wins
                    if row1[COL_IEVALUE] < row2[COL_IEVALUE]:
                        to_drop.add(idx2)
                    elif row2[COL_IEVALUE] < row1[COL_IEVALUE]:
                        to_drop.add(idx1)
                    else:
                        # tie — higher coverage wins
                        if row1[COL_COVERAGE] >= row2[COL_COVERAGE]:
                            to_drop.add(idx2)
                        else:
                            to_drop.add(idx1)

        kept.append(group[~group.index.isin(to_drop)])

    return pd.concat(kept).reset_index(drop=True)


def main(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", header=None)

    # convert numeric columns
    for col in [COL_IEVALUE, COL_ALI_FROM, COL_ALI_TO, COL_COVERAGE]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    print(f"Total rows before: {len(df)}")

    df_resolved = resolve_overlaps(df, overlap_threshold=0.30)

    print(f"Total rows after:  {len(df_resolved)}")
    print(f"Rows removed:      {len(df) - len(df_resolved)}")

    df_resolved.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python resolve_overlaps.py input.tsv output.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

######## sanity check to see how many were dropped ##############
for enzyme in GH29 GH95 GH141 GH107 GH168; do
    count1=$(awk -v e="$enzyme" '$1 ~ e {print $4}' prokaryotic_MAGs_vs_dbCAN.domtblout.tsv | sort -u | wc -l)
    count2=$(awk -v e="$enzyme" '$1 ~ e {print $4}' coverage_60_percent_domtblout.tsv | sort -u | wc -l)
    count3=$(awk -v e="$enzyme" '$1 ~ e {print $4}' domain_resolved_coverage_domtblout.tsv | sort -u | wc -l)
    echo "$enzyme  beforeanyfiltering: $count1  aftercoverage: $count2 coverage&overlapcleaned: $count3"

'''
GH29  beforeanyfiltering: 859  aftercoverage: 789 coverage&overlapcleaned: 789
GH95  beforeanyfiltering: 476  aftercoverage: 396 coverage&overlapcleaned: 390
GH141  beforeanyfiltering: 307  aftercoverage: 247 coverage&overlapcleaned: 247
GH107  beforeanyfiltering: 41  aftercoverage: 32 coverage&overlapcleaned: 32
GH168  beforeanyfiltering: 131  aftercoverage: 123 coverage&overlapcleaned: 123
'''







#################################### Protein Structure predictions ##############################
# esmfold basic prediction 
singularity exec \
  -B /scratch/shrinivas \
  --pwd /scratch/shrinivas/Sargassum/Degradation_experiment/metagenomic_analysis/substrate_predictions/cazyme_witch_hunt/gh29_esmfold \
  /scratch/singularity/ESMFold_v2.0.1_cuda_v11.8.0-rev1.sif \
  esm-fold \
  --fasta gh29_sequences.faa \
  -o pdb_output_files \
  -m /model \
  --cpu-only




