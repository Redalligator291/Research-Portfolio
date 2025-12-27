#!/bin/bash

# Specify the output file name
output_file="protein.ndx"

# Start and end numbers
start_number=1881
end_number=4666

# Create the index file
echo "[ Protein ]" > $output_file

# Loop through the numbers and add them to the index file
for ((i=start_number; i<=end_number; i++)); do
    echo "$i" >> $output_file
done

# Notify the user that the file has been created
echo "Protein index file created: $output_file"

