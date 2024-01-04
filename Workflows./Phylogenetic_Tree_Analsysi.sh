
# Mafft alignment

#!/bin/bash

input_dir="/input/directory"

# Set the output directory path
output_dir="enter/output/directory"
threads=40  # Set the number of threads

# Ensure the output directory exists
mkdir -p "$output_directory"

# Loop through each file in the input directory
for file in "$input_directory"/*.fa; do
    if [ -f "$file" ]; then
        # Extract the file name without extension
        filename=$(basename -- "$file")
        filename_no_extension="${filename%.*}"

        # Construct the output file name
        output_file="$output_directory/$filename_no_extension"_msa_align.fasta

        # Run the MAFFT command
        /home/shrinivas/Programs/mafft-linux64/mafft.bat --thread "$threads" --localpair --maxiterate 1000 --auto "$file" > "$output_file"

        echo "Alignment complete for $filename"
    fi
done


# Trimal loop
#!/bin/bash

# Set the input directory path
input_dir="/input/directory"

# Set the output directory path
output_dir="enter/output/directory "

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all files in the input directory with a .fasta extension
for file in "$input_dir"/*.fasta; do
    # Extract the filename without the extension
    filename=$(basename "$file" .fasta)

    # Define the output file path in the new directory
    output_file="$output_dir/${filename}_trim.fasta"

    # Run trimal command
    trimal -automated1 -in "$file" -out "$output_file"

    echo "Processed: $filename"
done


##### Iq tree #####

input_dir="/path/to/your/input/directory"

# Create a directory for the output files
output_dir="/path/to/your/output/directory"
mkdir -p "$output_dir"

for file in "$input_dir"/*.fasta; do
    # Extract the filename without the extension
    filename=$(basename "$file" .fasta)

    # Define the output file path in the new directory
    output_file="$output_dir/${filename}_iqtree"

    # Run iqtree command
    iqtree -s "$file" -m TEST -bb 1000 -nt AUTO -pre "$output_file"

    echo "Processed: $filename"
done
