#diamond_blast_orthogrups.sh

DIAMOND_PATH="/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond"

# Set the path to the query and database directories
QUERY_DIR="/scratch/shrinivas/Paulinella_Metagenome/Orthofinder_Analysis/Diamond_Blast_Orthogroups/Orthogroup_Sequences"
DB_DIR="/scratch/databases/ncbi/2022_07/nr.dmnd"

OUTPUT_DIR="/scratch/shrinivas/Paulinella_Metagenome/bin21"

# Set other parameters
EVALUE="0.001"
MAX_TARGET_SEQS="1"
THREADS= 50

# Loop through the files listed in files2run.txt
while read -r FILENAME; do
  # Generate the output filename based on the input filename
  QUERY="${QUERY_DIR}/${FILENAME}.faa"

  # Generate the output filename based on the input filename
  OUT="${OUTPUT_DIR}/${FILENAME}.blastp"

  # Run Diamond BLASTP
  ${DIAMOND_PATH} blastp --ultra-sensitive --max-target-seqs ${MAX_TARGET_SEQS} --evalue ${EVALUE} \
    --query "${QUERY_DIR}/${FILENAME}" --db "${DB_DIR}" --out "${OUT}" \
    --compress 1 --threads ${THREADS} --outfmt 6 qseqid sseqid pident length mismatch gapopen \
    qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms

done



##### for a single run ######
/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond blastp --ultra-sensitive --iterate --max-target-seqs 1 --evalue 0.001 --query /scratch/shrinivas/Paulinella_Metagenome/Oct24_2023/Outgroup_included_analysis/OrthoFinder/Results_Oct24_1/extracted_sequences.fa --db /scratch/databases/ncbi/2022_07/nr.dmnd --out Orthogroup_representative_BLAST.tsv --threads 70 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms


 /home/timothy/programs/DIAMOND_v2.1.8/bin/diamond blastx --ultra-sensitive --max-target-seqs 5 --evalue 0.00001 --query extracted_sequences.fasta --db /scratch/databases/ncbi/2022_07/nr.dmnd --out Extracted_sequences_non_paulinella_cadidates_BLAST.tsv --threads 20 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms



/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond blastp --ultra-sensitive --iterate --max-target-seqs 0 --evalue 0.00001 --query /scratch/shrinivas/Paulinella_Metagenome/Phylogenetic_Analysis_Orthogroups/blast_unlimited/extracted_sequences.fa --db /scratch/databases/ncbi/2022_07/nr.dmnd --out Phylogenetic_Analysis_Unlimited_BLAST.tsv --threads 70 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms



/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond blastp --ultra-sensitive --max-target-seqs 1 --evalue 0.00001 --query concatenated_sequences.fa --db /scratch/databases/ncbi/2022_07/nr.dmnd --out upset_orthogroup_analysis_blast_result.tsv --threads 60 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms
