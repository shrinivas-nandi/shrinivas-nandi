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

grep -v '^#' allMAGs_vs_dbCAN.domtblout | awk '
{
  print $4, $1, $13, $14, $7, $18, $19, $16, $17
}' OFS="\t" > domtblout.parsed.tsv


#!/usr/bin/env python3
import sys
from collections import defaultdict

OVERLAP_FRAC = 0.30
MIN_COVERAGE = 0.35

# overlap function: look for coordinates. are domains on the same part of the protein i.e., is this an hmm clash
def overlap(a_start, a_end, b_start, b_end):
    ov = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    a_len = a_end - a_start + 1
    b_len = b_end - b_start + 1
    return ov / min(a_len, b_len) if min(a_len, b_len) > 0 else 0

proteins = defaultdict(list)

with open(sys.argv[1]) as f:
    for line in f:
        (
            protein, hmm, ie, score, seqlen,
            s_start, s_end, h_start, h_end
        ) = line.strip().split("\t")

        hmm_len = abs(int(h_end) - int(h_start)) + 1 # generate hmm coverage feature
        coverage = hmm_len / int(seqlen) # calculate the coverage

        if coverage < MIN_COVERAGE:
            continue  # drop weak partial matches early

        proteins[protein].append({
            "hmm": hmm,
            "ie": float(ie),
            "score": float(score),
            "start": int(s_start),
            "end": int(s_end),
            "coverage": coverage
        })

def resolve(hits):
    kept = []
    used = set()

    for i, h1 in enumerate(hits):
        if i in used:
            continue

        group = [i]
        for j, h2 in enumerate(hits):
            if i != j:
                if overlap(h1["start"], h1["end"], h2["start"], h2["end"]) >= OVERLAP_FRAC:
                    group.append(j)

        for g in group:
            used.add(g)

        best = min(
            (hits[g] for g in group),
            key=lambda x: (x["ie"], -x["score"])
        )
        kept.append(best)

    return kept

print("protein\thmm\ti_evalue\tscore\tcoverage\tstart\tend")

for protein, hits in proteins.items():
    for h in resolve(hits):
        print(
            protein,
            h["hmm"],
            h["ie"],
            h["score"],
            f"{h['coverage']:.2f}",
            h["start"],
            h["end"],
            sep="\t"
        )

python resolve_domains_with_coverage.py domtblout.parsed.tsv > resolved_domains.tsv


#### select gh29 families###
# cut offs: coverage >= 0.35 (need to calculate on my own), 
awk '
$1 == "GH29" && $0 !~ /^#/ {
  hmm_len = $3
  cov = ($17 - $16 + 1) / hmm_len
  if (cov >= 0.35) print $4
}' allMAGs_vs_dbCAN.domtblout | sort -u > GH29_proteins.txt




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




