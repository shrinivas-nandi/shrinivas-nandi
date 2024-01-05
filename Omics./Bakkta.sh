for fasta_gz_file in *.fasta.gz; do
    # Check if the file exists before proceeding
    if [ -e "$fasta_gz_file" ]; then
        # Run the bakta tool on the gzipped file
        bakta --db /scratch/timothy/databases/bakta_db_v4.0/db --threads 10 "$fasta_gz_file"
    else
        echo "File not found: $fasta_gz_file"
    fi
done
