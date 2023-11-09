BLAST_nr.sh
# path to files /scratch/shrinivas/Paulinella_Metagenome/Prodigal/Protein_Files

DIAMOND_PATH="/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond"

# Specify the input protein FASTA file
INPUT_FASTA="Paulinella_Genome_Combined.faa"

# Specify the name for the DIAMOND database
DIAMOND_DB_NAME="Paulinella_Genome_Combined.fasta_DB"

# Run DIAMOND makedb to create the database
/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond makedb --in "Paulinella_Genome_Combined.faa" -d "Paulinella_Genome_Combined.fasta_DB" --dbtype prot





#Need to diamond blast against nr 

makeblastdb -in Paulinella_Genome_Combined.faa  -dbtype prot -out Paulinella_Genome_Combined.fasta_DB

#Go for top hit 

DIAMOND_PATH="/home/timothy/programs/DIAMOND_v2.1.8/bin/diamond"

# Set the path to the query and database directories
QUERY_DIR="/scratch/shrinivas/Paulinella_Metagenome/Prodigal/Protein_Files"
DB_DIR="/scratch/databases/ncbi/2022_07/nr.dmnd"
OUTPUT_DIR="/scratch/shrinivas/Paulinella_Metagenome/Prodigal/Protein_Files/BLAST_nr_output"

# Set other parameters
EVALUE="0.001"
MAX_TARGET_SEQS="1"
THREADS=50  # Removed the extra space after the equal sign

# Use a for loop to iterate over all .faa files in the directory
for FILE in "$QUERY_DIR"/*.faa; do
  # Extract the filename without the path and extension
  FILENAME=$(basename "$FILE" .faa)

  # Generate the output filename based on the input filename
  OUT="${OUTPUT_DIR}/${FILENAME}.blastp"

  # Run Diamond BLASTP
  "${DIAMOND_PATH}" blastp --ultra-sensitive --max-target-seqs "${MAX_TARGET_SEQS}" --evalue "${EVALUE}" \
    --query "$FILE" --db "${DB_DIR}" --out "$OUT" \
    --compress 1 --threads "${THREADS}" --outfmt 6 qseqid sseqid pident length mismatch gapopen \
    qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames sphylums skingdoms sskingdoms
done
