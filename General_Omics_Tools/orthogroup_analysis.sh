directory="/scratch/shrinivas/Paulinella/Paulinella_Metagenome/Paulinella_Genomes/Proteome/ovalis_only_genes/Potential_HGT/fasta_euk_comp_HGT/Orthogroup_Sequences"

# List file containing sequence names
list_file="gene_list.txt"

# Output TSV file
output_file="output.tsv"

# Clear existing output file
> "$output_file"

# Loop through each sequence name in the list file
while IFS= read -r sequence_name; do
            # Loop through each fasta file in the directory
                for file in "$directory"/*.fa; do
                                # Count the number of sequences in the file (excluding headers)
                                        sequence_count=$(grep -c "^>" "$file")

                                                # Check if the file contains at least two sequences
                                                        if [ "$sequence_count" -ge 2 ]; then
                                                                            # Check if the file contains the specified sequence header
                                                                                        if grep -q "^>$sequence_name" "$file"; then
                                                                                                                # Output filename and sequence header to TSV file
                                                                                                                                echo -e "$(basename "$file")\t$sequence_name" >> "$output_file"
                                                                                                                                                break
                                                                                                                                                            fi

     fi

         done

 done < "$list_file"
~                               
