#!/bin/bash

#SBATCH --output /nobackup/lab_menche/mmeyenberg/apobex/log_%j.log
#SBATCH --job-name=apobex
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load BEDTools/2.30.0-GCC-11.3.0

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 positions.txt"
    exit 1
fi

positions_file="$1"
reference_genome= </ref_genome/ref_genome.fasta>

# Ensure the positions file exists
if [ ! -f "$positions_file" ]; then
    echo "Error: Positions file '$positions_file' not found."
    exit 1
fi

# Ensure the reference genome file exists
if [ ! -f "$reference_genome" ]; then
    echo "Error: Reference genome file '$reference_genome' not found."
    exit 1
fi

# Generate the output file name
base_name=$(basename "$positions_file" .txt)
output_file="${base_name}_sequences.txt"

# Convert positions to BED format for bedtools
awk -v OFS="\t" '{ print $1, $2-21, $2+20 }' "$positions_file" > "${base_name}_positions.bed"

# Extract sequences using bedtools
bedtools getfasta -fi "$reference_genome" -bed "${base_name}_positions.bed" -fo "$output_file" -name

# Clean up intermediate BED file
rm "${base_name}_positions.bed"