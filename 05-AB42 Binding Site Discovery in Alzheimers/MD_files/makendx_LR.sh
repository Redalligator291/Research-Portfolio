#!/bin/bash

# Specify the output file name
output_file="protein.ndx"

# Start and end numbers for Receptor
start_receptor=1882
end_receptor=1952

# Create the index file
echo "[ Ligand ]" > $output_file

# Loop through the numbers and add them to the index file for Receptor
for ((i=start_receptor; i<=end_receptor; i++)); do
    echo "$i" >> $output_file
done


# Start and end numbers for Ligand
start_ligand=1
end_ligand=1881

# Append to the index file
echo "[ Receptor ]" >> $output_file

# Loop through the numbers and add them to the index file for Ligand
for ((i=start_ligand; i<=end_ligand; i++)); do
    echo "$i" >> $output_file
done

# Notify the user that the file has been created
echo "Protein index file created: $output_file"

