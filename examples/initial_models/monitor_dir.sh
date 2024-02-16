#!/bin/bash

# Check if input folder argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_folder>"
    exit 1
fi

input_folder="$1"

# Check if input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Error: Input folder '$input_folder' does not exist."
    exit 1
fi

# Loop through files in input folder
for file in "$input_folder"/*; do
    # Check if file is a regular file
    if [ -f "$file" ]; then
        # Get filename without extension
        filename=$(basename -- "$file")
        filename_no_ext="${filename%.*}"
        
        # Execute command with nohup
        nohup julia monitor.jl "$filename_no_ext" &
        
        echo "Started monitoring $filename_no_ext"
    fi
done

echo "Monitoring started for all files in $input_folder"