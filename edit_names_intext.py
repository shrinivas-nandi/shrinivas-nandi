# this is a python script to edit names in a file 
# like a tree file or something like that 

def replace_ep_values(tree_filename, replacements_filename, output_filename):
    """Reads the input tree file and the list of replacements, performs the substitution, and writes the output."""
    # Read the tree file
    with open(tree_filename, 'r') as tree_file:
        tree_content = tree_file.read()

    # Read the replacement list
    with open(replacements_filename, 'r') as replacements_file:
        replacement_lines = replacements_file.readlines()

    # Create a mapping of EP values to replacement values
    replacement_mapping = {}
    for line in replacement_lines:
        line = line.strip().split('_')
        ep_value = line[0]
        replacement_value = line[1].split('.')[0]
        replacement_mapping[ep_value] = replacement_value

    # Perform the substitution
    for ep_value, replacement_value in replacement_mapping.items():
        tree_content = tree_content.replace(ep_value, replacement_value)

    # Write the modified tree to the output file
    with open(output_filename, 'w') as output_file:
        output_file.write(tree_content)


# Replace these filenames with your actual filenames
tree_file_path = 'your_tree_file.txt'
replacements_file_path = 'your_replacements_file.txt'
output_file_path = 'output_tree.txt'

# Call the function
replace_ep_values(tree_file_path, replacements_file_path, output_file_path)
