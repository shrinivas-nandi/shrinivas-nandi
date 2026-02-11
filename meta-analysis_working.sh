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
/home/timothy/programs/iqtree-1.6.12-Linux/bin/iqtree -s /scratch/shrinivas/Sargassum/meta-analysis/naive_atlas_output/genomes/prokayotic_genomes/gtdbtk_identify/mafft_output/align/gtdbtk.bac120.user_msa.aligned.trimmed.fa -m TEST -bb 1000 -nt AUTO -pre /scratch/shrinivas/Sargassum/meta-analysis/naive_atlas_output/genomes/prokayotic_genomes/gtdbtk_identify/mafft_output/align/iqtree_results

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

1,581 highly likely cazymes in 284 mags from 1338 ~ 21% of MAGs
'''

###### repeat for alginate lyases #######
for enzyme in PL5 PL6 PL7 PL8 PL14 PL15 PL17 PL18 PL31 PL32 PL34 PL36 PL39 PL41; do
    count1=$(awk -v e="$enzyme" '$1 ~ e {print $4}' prokaryotic_MAGs_vs_dbCAN.domtblout.tsv | sort -u | wc -l)
    count2=$(awk -v e="$enzyme" '$1 ~ e {print $4}' coverage_60_percent_domtblout.tsv | sort -u | wc -l)
    count3=$(awk -v e="$enzyme" '$1 ~ e {print $4}' use_this_domain_resolved_coverage_domtblout.tsv | sort -u | wc -l)
    echo "$enzyme  beforeanyfiltering: $count1  aftercoverage: $count2 coverage&overlapcleaned: $count3"
done

'''
PL5  beforeanyfiltering: 26  aftercoverage: 24 coverage&overlapcleaned: 23
PL6  beforeanyfiltering: 730  aftercoverage: 670 coverage&overlapcleaned: 665
PL7  beforeanyfiltering: 1089  aftercoverage: 1048 coverage&overlapcleaned: 1048
PL8  beforeanyfiltering: 41  aftercoverage: 40 coverage&overlapcleaned: 39
PL14  beforeanyfiltering: 207  aftercoverage: 202 coverage&overlapcleaned: 198
PL15  beforeanyfiltering: 143  aftercoverage: 142 coverage&overlapcleaned: 141
PL17  beforeanyfiltering: 415  aftercoverage: 413 coverage&overlapcleaned: 399
PL18  beforeanyfiltering: 47  aftercoverage: 47 coverage&overlapcleaned: 47
PL31  beforeanyfiltering: 88  aftercoverage: 83 coverage&overlapcleaned: 81
PL32  beforeanyfiltering: 0  aftercoverage: 0 coverage&overlapcleaned: 0
PL34  beforeanyfiltering: 33  aftercoverage: 32 coverage&overlapcleaned: 32
PL36  beforeanyfiltering: 168  aftercoverage: 168 coverage&overlapcleaned: 4
PL39  beforeanyfiltering: 135  aftercoverage: 78 coverage&overlapcleaned: 13
PL41  beforeanyfiltering: 43  aftercoverage: 38 coverage&overlapcleaned: 2

3034 unique alginate lyases 
distributed in over 730 pMAGs (~ 54%)
'''
#################################### pull out the whole protein sequence ##############################
# extract proteins based on protein name 
awk '{print $4}' fucoidanases.tsv | sort -u > fucoidanase_protein_names.txt

/home/timothy/programs/seqkit_v2.3.1/seqkit grep -f fucoidanase_protein_names.txt /scratch/shrinivas/Sargassum/meta-analysis/predicted_proteins/prokaryotic_proteome/prokaryotic_proteome.faa > fucoidanases_full_protein.fasta

# repeat for alginates
awk '{print $4}' alginates.tsv | sort -u > alginate_protein_names.txt

/home/timothy/programs/seqkit_v2.3.1/seqkit grep -f alginate_protein_names.txt /scratch/shrinivas/Sargassum/meta-analysis/predicted_proteins/prokaryotic_proteome/prokaryotic_proteome.faa > alginate_lyases_full_protein.fasta

# extract catalytic activity using the ali to ali for each protein 
## we will use seqio to extract just the catalytic domain 

from Bio import SeqIO
import pandas as pd

df = pd.read_csv("/scratch/shrinivas/Sargassum/meta-analysis/predicted_proteins/prokaryotic_proteome/cazyme_prokaryotes/fucoidanases/fucoidanases.tsv", sep="\t", header=None)
seqs = SeqIO.index("/scratch/shrinivas/Sargassum/meta-analysis/predicted_proteins/prokaryotic_proteome/cazyme_prokaryotes/fucoidanases/fucoidanases_full_protein.fasta, "fasta")

with open("fucoidanases_gh_domain_only.fasta", "w") as out:
    for _, row in df.iterrows():
        protein = row[3]
        ali_from = int(row[17]) - 1  # convert to 0-indexed
        ali_to = int(row[18])
        if protein in seqs:
            subseq = seqs[protein].seq[ali_from:ali_to]
            out.write(f">{protein}_{ali_from+1}_{ali_to}\n{subseq}\n")


### verify with another tool, blast maybe#####

### pfamscan the rest for other domains
/scratch/shrinivas/programs/my_interproscan/interproscan-5.77-108.0/interproscan.sh \
  -i prokaryotic_plasmid_full_proteins.fasta \
  -f tsv \
  -b pfam_annotation_fullproteins \
  -dp \
  -appl Pfam \
  --cpu 60


####################### PLASMID PROTEINS ##################################
hmmscan \
  --cpu 64 \
  --domtblout plasmid_MAGs_vs_dbCAN.domtblout \
  --tblout plasmid_MAGs_vs_dbCAN.tblout \
  -E 1e-15 --domE 1e-15 \
  /scratch/shrinivas/Databases/CaZy_db/dbCAN-HMMdb-V14.txt \
  /scratch/shrinivas/Sargassum/meta-analysis/predicted_proteins/prokaryotic_proteome/plasmid_proteins.faa \
  > all_PLASMID_MAGs_vs_dbCAN.log


for enzyme in PL5 PL6 PL7 PL8 PL14 PL15 PL17 PL18 PL31 PL32 PL34 PL36 PL39 PL41; do
    count1=$(awk -v e="$enzyme" '$1 ~ e {print $4}' plasmid_MAGs_vs_dbCAN.domtblout | sort -u | wc -l)
done

# also completed for fucoidan (only 1 GH29)
oveall very little stuff in the plasmids. 

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




