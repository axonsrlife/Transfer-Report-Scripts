#!/bin/bash

# Change letters to indicate charges
change_to_charge_indicator() {
  sequence=$1
  result=""

  for ((i=0; i<${#sequence}; i++)); do
    amino_acid="${sequence:$i:1}"

    case $amino_acid in
      "A"|"G"|"V"|"C"|"P"|"L"|"I"|"M"|"W"|"F"|"S"|"T"|"Y"|"N"|"Q") charge=0;;
      "K"|"R"|"H") charge=1;;
      "D"|"E") charge=-1;;
      *) charge="unknown";;
    esac

    result="${result}${charge} "
  done

  echo "$result"
}

# Check if arguments are provided
if [ $# -lt 2 ]; then
  echo "Usage: $0 <input_file> <output_file>"
  exit 1
fi

# Get the input and output files from the command line arguments
input_file=$1
output_file=$2

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 1
fi

# Process each sequence in the input file and write the results to the output file
while read -r line; do
  # If the line starts with ">", it's a header line
  if [[ $line == ">"* ]]; then
    # Print the header to the output file
    echo "$line" >> "$output_file"
  else
    # Process the sequence and write the charge indicator to the output file
    result=$(change_to_charge_indicator "$line")
    echo "$result" >> "$output_file"
  fi
done < "$input_file"
