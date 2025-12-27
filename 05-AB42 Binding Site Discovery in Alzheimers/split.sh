#!/bin/bash

# Input file name
input_file="result.dock4"

# Initialize variables
output_file=""
cluster_num=""
cluster_member=""

# Read the input file line by line
while IFS= read -r line; do
    # Check for REMARK lines to extract CLUSTER_NUM and CLUSTER_MEMBER
    if [[ $line == "REMARK CLUSTER_NUM "* ]]; then
        cluster_num=$(echo "$line" | awk '{print $4}')
    elif [[ $line == "REMARK CLUSTER_MEMBER "* ]]; then
        cluster_member=$(echo "$line" | awk '{print $4}')
        output_file="${cluster_num}_${cluster_member}.pdb"
    fi
    
    # Check for TER line to close the current file
    if [[ $line == "TER" ]]; then
        echo "$line" >> "$output_file"
        output_file="" # Reset output file name
        continue
    fi
    
    # Write the line to the current output file if it's defined
    if [[ -n $output_file ]]; then
        echo "$line" >> "$output_file"
    fi

done < "$input_file"

echo "Files created successfully."
