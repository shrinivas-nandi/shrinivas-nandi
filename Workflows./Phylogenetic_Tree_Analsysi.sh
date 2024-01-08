
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



####### 
#######
# This bit of code here is done on R. 
# Effectively a loop for plotting out all the tree files. 
library(ape)

# Function to assign colors based on species names
assign_tip_color <- function(species_name) {
  if (grepl("KR01", species_name)) {
    return("blue")
  } else if (grepl("^g", species_name)) {
    return("red")
  } else if (grepl("CCAC0185", species_name)) {
    return("green")
  } else {
    return("orange")  # Default color for other species
  }
}

# Function to plot and save tree as PDF
plot_and_save_tree <- function(tree_file_path, output_pdf_path) {
  # Read the tree
  tree <- read.tree(tree_file_path)
  
  # Assign tip colors
  tip_colors <- sapply(tree$tip.label, assign_tip_color)
  
  # Plot the tree with adjusted width and height
  pdf(output_pdf_path, width = 20, height = 12)  # Adjust the width and height as needed
  plot(tree, edge.width = 1.5, cex = 1, tip.color = tip_colors)
  nodelabels(tree$node.label, node = 2:(tree$Nnode + Ntip(tree)), adj = c(1, -0.2), frame = "none", cex = 0.5)
  dev.off()
}

# Directory containing tree files
tree_directory <- "/Users/shrinivas/Desktop/Paulinella/P_ovalis_Metagenome/Tree_plot/"

# List all tree files in the directory
tree_files <- list.files(tree_directory, pattern = "\\.treefile$", full.names = TRUE)

# Loop through each tree file
for (tree_file in tree_files) {
  # Generate output PDF path based on the input tree file name
  output_pdf_path <- paste0(tools::file_path_sans_ext(tree_file), "_colored.pdf")
  
  # Plot and save the tree with adjusted width and height
  plot_and_save_tree(tree_file, output_pdf_path)
}


